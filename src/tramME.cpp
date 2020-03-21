/*==================================================================
  Transformation model: F_Y(y|x,gamma) = F_Z(h(y|x,gamma))
  Structure of ME: h(y|x,gamma) = (a(y)', b(x)') * beta + Z * gamma
  ==================================================================
  TODO:
  - truncation
  - add more complex random effects
 */

#include <TMB.hpp>

// Model types
enum valid_modtype {
    Lm = 0
};

// Valid error distributions
enum valid_errdist {
    Normal = 0, Logistic = 1, MinExtrVal = 2, MaxExtrVal = 3
};

// F_Z(h(y)), cens: censoring indicator, if not finite, return value
template <class Type>
Type cdf(Type hy, int errdist) {
  Type out;
  switch (errdist) {
    case Normal:
      out = pnorm(hy);
      break;
    case Logistic:
      //out = Type(1) / (Type(1) + exp(-hy));
      out = invlogit(hy);
      break;
    case MinExtrVal:
      out = Type(1) - exp(-exp(hy));
      break;
    case MaxExtrVal:
      out = exp(-exp(-hy));
      break;
    default:
      error("Unknown error distribution!");
  }
  return out;
}

// log-density error function
template <class Type>
Type ldens(Type hy, int errdist) {
  Type out;
  switch (errdist) {
    case Normal:
      out = dnorm(hy, Type(0), Type(1), true);
      break;
    case Logistic:
      out = dlogis(hy, Type(0), Type(1), true);
      break;
    case MinExtrVal:
      out = hy - exp(hy);
      break;
    case MaxExtrVal:
      out = - hy - exp(-hy);
      break;
    default:
      error("Unknown error distribution!");
  }
  return out;
}

GVECTORIZE(ldens, V, I, N, N, N, N);

// Covariance terms of the random effects
template <class Type>
struct re_cov_term {
  vector<Type> sd;
  matrix<Type> corr;
};

// Negative log-density of the random effects
template <class Type>
Type re_nldens(vector<Type> gamma, vector<Type> theta, int blocksize, re_cov_term<Type>& term) {
  Type ans = 0;
  if (blocksize == 1) { // diagonal cov matrix of the term
    //Type sd = theta[0];
    Type sd = exp(theta[0]); // parateterized w/ log(sd)
    ans -= dnorm(gamma, Type(0), sd, true).sum();
    //term.sd = theta;
    term.sd = exp(theta);
    matrix<Type> corr(1,1);
    term.corr = corr.setIdentity(1,1);
  } else { // correlated random effects (unstructured corr mat)
    int nblocks = gamma.size() / blocksize;
    //vector<Type> sd = theta.head(blocksize);
    vector<Type> sd = exp(theta.head(blocksize)); // parameterized w/ log(sd)
    vector<Type> corr_tr = theta.tail(theta.size() - blocksize);
    density::UNSTRUCTURED_CORR_t<Type> nldens(corr_tr);
    density::VECSCALE_t<density::UNSTRUCTURED_CORR_t<Type> > scnldens = density::VECSCALE(nldens, sd);
    for (int i = 0; i < nblocks; i++) {
      ans += scnldens(gamma.segment(i * blocksize, blocksize));
    }
    term.sd = sd;
    term.corr = nldens.cov();
  }
  return ans;
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_INTEGER(modtype); // model type for outputs
  DATA_INTEGER(errdist); // Type of error distribution
  // model matrices: left, right, interval (l & r), exact (also prime)
  DATA_MATRIX(MMl);
  DATA_MATRIX(MMr);
  DATA_MATRIX(MMil);
  DATA_MATRIX(MMir);
  DATA_MATRIX(MMe);
  DATA_MATRIX(MMeprime); // for constructing h(y)'
  // design matrices of the random effects
  DATA_SPARSE_MATRIX(Zl); // left-cens
  DATA_SPARSE_MATRIX(Zr); // right-cens
  DATA_SPARSE_MATRIX(Zi); // interval-cens
  DATA_SPARSE_MATRIX(Ze); // exact
  DATA_IVECTOR(re_termsize);  // number of REs corresponding to each term
  DATA_IVECTOR(re_blocksize); // size of blocks in the cov matrix of REs corresponding one term
  // offsets
  DATA_VECTOR(offsetl);
  DATA_VECTOR(offsetr);
  DATA_VECTOR(offseti);
  DATA_VECTOR(offsete);
  // weights
  DATA_VECTOR(weightsl);
  DATA_VECTOR(weightsr);
  DATA_VECTOR(weightsi);
  DATA_VECTOR(weightse);
  // parameters
  PARAMETER_VECTOR(beta);  // fixed effects
  PARAMETER_VECTOR(gamma); // random effects
  PARAMETER_VECTOR(theta); // covariance parameters

  parallel_accumulator<Type> nll(this);

  vector<re_cov_term<Type> > re_cov(re_termsize.size()); // for reporting purposes

  // ====== Likelihood contributions of RE terms
  int cum_ts = 0;
  int cum_bs = 0;
  for (int i = 0; i < re_termsize.size(); i++) {
    int nth = re_blocksize(i) * (re_blocksize(i)+1) / 2; // nr of covariance parameters
    vector<Type> gamma_s = gamma.segment(cum_ts, re_termsize(i));
    vector<Type> theta_s = theta.segment(cum_bs, nth);
    nll += re_nldens(gamma_s, theta_s, re_blocksize(i), re_cov(i));
    cum_ts += re_termsize(i);
    cum_bs += nth;
  }

  // ====== Likelihood contributions of observations
  // --- Left-censored observations
  if (MMl.rows() > 0) {
    vector<Type> hl = MMl * beta + Zl * gamma + offsetl;
    for (int i = 0; i < MMl.rows(); i++) {
      nll -= log(cdf(hl(i), errdist)) * weightsl(i);
    }
  }

  // --- Right-censored observations
  if (MMr.rows() > 0) {
    vector<Type> hr = MMr * beta + Zr * gamma + offsetr;
    for (int i = 0; i < MMr.rows(); i++) {
      nll -= log(Type(1) - cdf(hr(i), errdist)) * weightsr(i);
    }
  }

  // --- Interval-censored observations
  if (MMil.rows() > 0) {
    vector<Type> hil = MMil * beta + Zi * gamma + offseti;
    vector<Type> hir = MMir * beta + Zi * gamma + offseti;
    for (int i = 0; i < MMil.rows(); i++) {
      nll -= log(cdf(hir(i), errdist) - cdf(hil(i), errdist)) * weightsi(i);
    }
  }

  // --- Exact observations
  if (MMe.rows() > 0) {
    vector<Type> he = MMe * beta + Ze * gamma + offsete;
    vector<Type> hprime = MMeprime * beta;
    nll -= ((ldens(he, errdist) + log(hprime)) * weightse).sum();
  }

  // ====== Reporting SD and correlation matrices of random effects
  vector<matrix<Type> > corr_rep(re_cov.size());
  vector<vector<Type> > sd_rep(re_cov.size());
  for(int i = 0; i < re_cov.size(); i++) {
      corr_rep(i) = re_cov(i).corr;
      sd_rep(i) = re_cov(i).sd;
  }
  REPORT(corr_rep);
  REPORT(sd_rep);

  // ====== Report transformed paremeters for specific model types
  if (modtype == Lm) {
    int p = beta.size() - 1;
    Type sigma = Type(1) / beta(1);
    vector<Type> b(p);
    vector<Type> th = theta;
    // --- Coefficients
    b(0) = -beta(0) * sigma;
    b.tail(p-1) = beta.tail(p-1) * sigma;
    // --- Random effects parameters
    int cum_bs = 0;
    for (int i = 0; i < re_termsize.size(); i++) {
      int nth = re_blocksize(i) * (re_blocksize(i)+1) / 2; // nr of covariance parameters
      th.segment(cum_bs, re_blocksize(i)) += log(sigma);
      cum_bs += nth;
    }
    ADREPORT(b);
    ADREPORT(sigma);
    ADREPORT(th);
  }

  return nll;
}
