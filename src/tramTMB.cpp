/*==================================================================
  Transformation model: F_Y(y|x,gamma) = F_Z(h(y|x,gamma))
  Structure of ME: h(y|x,gamma) = (a(y)', b(x)') * beta + Z * gamma
  ==================================================================
 */

#include <TMB.hpp>


// Valid error distributions
enum valid_errdist {
    Normal = 0, Logistic = 1, MinExtrVal = 2, MaxExtrVal = 3, Exponential = 4
};

// Inverse link functions
template <class Type>
Type cdf(Type x, int errdist) {
  Type out;
  switch (errdist) {
    case Normal:
      out = pnorm(x);
      break;
    case Logistic:
      out = invlogit(x);
      break;
    case MinExtrVal:
      out = Type(1) - exp(-exp(x));
      break;
    case MaxExtrVal:
      out = exp(-exp(-x));
      break;
    case Exponential:
      out = pexp(x, Type(1));
      break;
    default:
      error("Unknown error distribution!");
  }
  return out;
}

// Log-density of the error distribution
template <class Type>
Type ldens(Type x, int errdist) {
  Type out;
  switch (errdist) {
    case Normal:
      out = dnorm(x, Type(0), Type(1), true);
      break;
    case Logistic:
      out = dlogis(x, Type(0), Type(1), true);
      break;
    case MinExtrVal:
      out = x - exp(x);
      break;
    case MaxExtrVal:
      out = - x - exp(-x);
      break;
    case Exponential:
      out = dexp(x, Type(1), true);
      break;
    default:
      error("Unknown error distribution!");
  }
  return out;
}

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
  DATA_IVECTOR(errdist); // Type of error distribution
  // model matrices: left, right, interval (l & r), exact (also prime)
  DATA_MATRIX(MMl);
  DATA_MATRIX(MMr);
  DATA_MATRIX(MMil);
  DATA_MATRIX(MMir);
  DATA_MATRIX(MMe);
  DATA_MATRIX(MMeprime); // for constructing h(y)'
  // indices
  DATA_IVECTOR(whichl);
  DATA_IVECTOR(whichr);
  DATA_IVECTOR(whichi);
  DATA_IVECTOR(whiche);
  DATA_SPARSE_MATRIX(Z); // design matrix of the random effects
  DATA_IVECTOR(re_termsize);  // number of REs corresponding to each term
  DATA_IVECTOR(re_blocksize); // size of blocks in the cov matrix of REs corresponding one term
  bool re_flag = (re_termsize.size() > 0);
  // offsets & weights (optionally updateable)
  DATA_VECTOR(offset);
  DATA_VECTOR(weights);
  DATA_INTEGER(do_update);
  bool upd_flag = (do_update > 0);
  if (upd_flag) {
    DATA_UPDATE(offset);
    DATA_UPDATE(weights);
  }
  // truncation
  DATA_MATRIX(MTl);
  DATA_MATRIX(MTr);
  DATA_MATRIX(MTil);
  DATA_MATRIX(MTir);
  DATA_IVECTOR(whichtl);
  DATA_IVECTOR(whichtr);
  DATA_IVECTOR(whichti);
  // parameters
  PARAMETER_VECTOR(beta);  // fixed effects
  PARAMETER_VECTOR(theta); // covariance parameters
  PARAMETER_VECTOR(alpha0); // auxiliary params for score residuals
  PARAMETER_VECTOR(gamma); // random effects
  bool resid_flag = (alpha0.size() > 0);

  parallel_accumulator<Type> nll(this);

  vector<re_cov_term<Type> > re_cov(re_termsize.size()); // for reporting purposes

  vector<Type> Zg(errdist.size()); // default: w/o random effects
  Zg.setZero();

  if (re_flag) {
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

    Zg = Z * gamma;
  }

  vector<Type> res(errdist.size()); // default: w/o auxiliary parameters
  res.setZero();

  if (resid_flag) {
    res += alpha0;
  }

  // ====== Likelihood contributions of observations
  // --- Left-censored observations
  if (whichl.size() > 0) {
    vector<Type> hl = MMl * beta + Zg(whichl) + offset(whichl) + res(whichl);
    for (int i = 0; i < whichl.size(); i++) {
      nll -= log(1.e-20 + cdf(hl(i), errdist(whichl(i)))) * weights(whichl(i));
    }
  }

  // --- Right-censored observations
  if (whichr.size() > 0) {
    vector<Type> hr = MMr * beta + Zg(whichr) + offset(whichr) + res(whichr);
    for (int i = 0; i < whichr.size(); i++) {
      nll -= log(1.e-20 + Type(1) - cdf(hr(i), errdist(whichr(i)))) * weights(whichr(i));
    }
  }

  // --- Interval-censored observations
  if (whichi.size() > 0) {
    vector<Type> hil = MMil * beta + Zg(whichi) + offset(whichi) + res(whichi);
    vector<Type> hir = MMir * beta + Zg(whichi) + offset(whichi) + res(whichi);
    for (int i = 0; i < whichi.size(); i++) {
      nll -= log(1.e-20 + cdf(hir(i), errdist(whichi(i))) - cdf(hil(i), errdist(whichi(i)))) *
        weights(whichi(i));
    }
  }

  // --- Exact observations
  if (whiche.size() > 0) {
    vector<Type> he = MMe * beta + Zg(whiche) + offset(whiche) + res(whiche);
    vector<Type> hprime = MMeprime * beta;
    for (int i = 0; i < whiche.size(); i++) {
      nll -= (ldens(he(i), errdist(whiche(i))) + log(hprime(i))) * weights(whiche(i));
    }
  }

  // --- Left truncation
  if (whichtl.size() > 0) {
    vector<Type> htl = MTl * beta + Zg(whichtl) + offset(whichtl) + res(whichtl);
    for (int i = 0; i < whichtl.size(); i++) {
      nll += log(1.e-20 + Type(1) - cdf(htl(i), errdist(whichtl(i)))) * weights(whichtl(i));
    }
  }

  // --- Right truncation
  if (whichtr.size() > 0) {
    vector<Type> htr = MTr * beta + Zg(whichtr) + offset(whichtr) + res(whichtr);
    for (int i = 0; i < whichtr.size(); i++) {
      nll += log(1.e-20 + cdf(htr(i), errdist(whichtr(i)))) * weights(whichtr(i));
    }
  }

  // --- Interval truncation
  if (whichti.size() > 0) {
    vector<Type> htil = MTil * beta + Zg(whichti) + offset(whichti) + res(whichti);
    vector<Type> htir = MTir * beta + Zg(whichti) + offset(whichti) + res(whichti);
    for (int i = 0; i < whichti.size(); i++) {
      nll += log(1.e-20 + cdf(htir(i), errdist(whichti(i))) - cdf(htil(i), errdist(whichti(i)))) *
        weights(whichti(i));
    }
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

  return nll;
}
