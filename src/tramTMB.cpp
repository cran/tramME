/*==================================================================
  Mixed-effects linear transformation models
  ------------------------------------------
  - various (at least ordered) outcomes
  - outcome error distributions may differ on the observation level
    (to be used with joint models)
  - exactly observed, randomly censored (left, right, interval) or
    truncated observations
  - the different inverse link functions imply different shift
    coef interpretations
  - random effects are optional
  - residuals calculated as the gradients of an additional intercept
    parameter (alpha0) fixed at 0
  - offsets and weights are optionally updateable for efficient
    calculations in certain potential extensions (e.g. trees and
    model-based boosting)
    TODO: post-estimation calculations on various scales
  ==================================================================
 */

#include <TMB.hpp>


// Valid error distributions
enum valid_errdist {
    Normal = 0, Logistic = 1, MinExtrVal = 2, MaxExtrVal = 3, Exponential = 4
};

enum valid_postproc {
  lp = 1, trafo = 2
};

enum valid_part {
  post = 0, lik = 1, pri = 2
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
  // === Inverse link functions (might be different for individual observations)
  DATA_IVECTOR(errdist); // Type of error distribution
  // === Baseline & interacting terms
  // -- Censoring
  // design matrices: left, right, interval (l & r), exact (also prime)
  DATA_MATRIX(Yl);
  DATA_MATRIX(Yr);
  DATA_MATRIX(Yil);
  DATA_MATRIX(Yir);
  DATA_MATRIX(Ye);
  DATA_MATRIX(Yeprime); // for constructing h(y)'
  // indices
  DATA_IVECTOR(whichl);
  DATA_IVECTOR(whichr);
  DATA_IVECTOR(whichi);
  DATA_IVECTOR(whiche);
  // -- Truncation
  DATA_MATRIX(Ytl);
  DATA_MATRIX(Ytr);
  DATA_MATRIX(Ytil);
  DATA_MATRIX(Ytir);
  // indices
  DATA_IVECTOR(whichtl);
  DATA_IVECTOR(whichtr);
  DATA_IVECTOR(whichti);
  // -- Coefficients
  PARAMETER_VECTOR(beta0);  // coef of baseline & interacting
  // === Fixed effects shift terms
  DATA_MATRIX(X);
  PARAMETER_VECTOR(beta);  // fixed effects
  // === Random effects
  DATA_SPARSE_MATRIX(Z); // design matrix of the random effects
  DATA_IVECTOR(re_termsize);  // number of REs corresponding to each term
  DATA_IVECTOR(re_blocksize); // size of blocks in the cov matrix of REs corresponding one term
  PARAMETER_VECTOR(gamma); // random effects vector
  PARAMETER_VECTOR(theta); // covariance parameters
  // === Offsets & weights (optionally updateable)
  DATA_VECTOR(offset);
  DATA_VECTOR(weights);
  DATA_INTEGER(do_update);
  if (do_update) {
    DATA_UPDATE(offset);
    DATA_UPDATE(weights);
  }
  // === Auxiliary params for score residuals
  PARAMETER_VECTOR(alpha0);
  // === Post-estimation calculations
  DATA_INTEGER(as_lm);
  DATA_INTEGER(postest_scale);
  DATA_MATRIX(Ype)
  DATA_MATRIX(Xpe);
  DATA_SPARSE_MATRIX(Zpe);
  // FIXME: is offset needed? DATA_VECTOR(offsetpe);
  // === Evaluate different parts of the negative ll
  DATA_INTEGER(part); // 0 = posterior, 1 = likelihood, 2 = prior

  Type npld = 0; // negative prior log-density
  parallel_accumulator<Type> nll(this);  // negative log-likelihood


  // ====== Fixed effects shift terms
  vector<Type> Xb(errdist.size()); // default: w/o random effects
  Xb.setZero();
  if (X.cols()) {
    Xb = X * beta;
  }

  // ====== Likelihood contributions of RE terms
  vector<re_cov_term<Type> > re_cov(re_termsize.size()); // for reporting purposes

  vector<Type> Zg(errdist.size()); // default: w/o random effects
  Zg.setZero();

  // if (re_flag) {
  if (re_termsize.size()) {
    int cum_ts = 0;
    int cum_bs = 0;
    for (int i = 0; i < re_termsize.size(); i++) {
      int nth = re_blocksize(i) * (re_blocksize(i)+1) / 2; // nr of covariance parameters
      vector<Type> gamma_s = gamma.segment(cum_ts, re_termsize(i));
      vector<Type> theta_s = theta.segment(cum_bs, nth);
      npld += re_nldens(gamma_s, theta_s, re_blocksize(i), re_cov(i));
      cum_ts += re_termsize(i);
      cum_bs += nth;
    }

    Zg = Z * gamma;
  }

  // ====== Residual calculations
  vector<Type> res(errdist.size()); // default: w/o auxiliary parameters
  res.setZero();

  // if (resid_flag) {
  if (alpha0.size()) {
    res += alpha0;
  }

  // ====== Likelihood contributions of observations
  // --- Left-censored observations
  if (whichl.size()) {
    vector<Type> hl = Yl * beta0 + Xb(whichl) + Zg(whichl) + offset(whichl) + res(whichl);
    for (int i = 0; i < whichl.size(); i++) {
      nll -= log(1.e-20 + cdf(hl(i), errdist(whichl(i)))) * weights(whichl(i));
    }
  }

  // --- Right-censored observations
  if (whichr.size()) {
    vector<Type> hr = Yr * beta0 + Xb(whichr) + Zg(whichr) + offset(whichr) + res(whichr);
    for (int i = 0; i < whichr.size(); i++) {
      nll -= log(1.e-20 + Type(1) - cdf(hr(i), errdist(whichr(i)))) * weights(whichr(i));
    }
  }

  // --- Interval-censored observations
  if (whichi.size()) {
    vector<Type> si(whichi.size());
    si = Xb(whichi) + Zg(whichi) + offset(whichi) + res(whichi);
    vector<Type> hil = Yil * beta0 + si;
    vector<Type> hir = Yir * beta0 + si;
    for (int i = 0; i < whichi.size(); i++) {
      nll -= log(1.e-20 + cdf(hir(i), errdist(whichi(i))) - cdf(hil(i), errdist(whichi(i)))) *
        weights(whichi(i));
    }
  }

  // --- Exact observations
  if (whiche.size()) {
    vector<Type> he = Ye * beta0 + Xb(whiche) + Zg(whiche) + offset(whiche) + res(whiche);
    vector<Type> hprime = Yeprime * beta0;
    for (int i = 0; i < whiche.size(); i++) {
      nll -= (ldens(he(i), errdist(whiche(i))) + log(hprime(i))) * weights(whiche(i));
    }
  }

  // --- Left truncation
  if (whichtl.size()) {
    vector<Type> htl = Ytl * beta0 + Xb(whichtl) + Zg(whichtl) + offset(whichtl) + res(whichtl);
    for (int i = 0; i < whichtl.size(); i++) {
      nll += log(1.e-20 + Type(1) - cdf(htl(i), errdist(whichtl(i)))) * weights(whichtl(i));
    }
  }

  // --- Right truncation
  if (whichtr.size()) {
    vector<Type> htr = Ytr * beta0 + Xb(whichtr) + Zg(whichtr) + offset(whichtr) + res(whichtr);
    for (int i = 0; i < whichtr.size(); i++) {
      nll += log(1.e-20 + cdf(htr(i), errdist(whichtr(i)))) * weights(whichtr(i));
    }
  }

  // --- Interval truncation
  if (whichti.size()) {
    vector<Type> sti(whichti.size());
    sti = Xb(whichti) + Zg(whichti) + offset(whichti) + res(whichti);
    vector<Type> htil = Ytil * beta0 + sti;
    vector<Type> htir = Ytir * beta0 + sti;
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

  // ====== Post-estimation calculations
  // --- Reparameterize as an LMM
  if (as_lm) {
    int p = 1 + beta.size();
    Type sigma = Type(1) / beta0(1);
    vector<Type> b(p);
    vector<Type> th = theta;
    // --- Coefficients
    b(0) = -beta0(0) * sigma;
    b.tail(p-1) = beta * sigma;
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

  if (postest_scale > 0) {
    // TODO: ADREPORT post-estimation calculations
    switch (postest_scale) {
      case lp: {
        vector<Type> pred = Xpe * beta + Zpe * gamma;
        if (as_lm) {
          Type sigma = Type(1) / beta0(1);
          pred *= sigma;
        }
        ADREPORT(pred);
        } break;
      case trafo: {
        vector<Type> pred = Ype * beta0 + Xpe * beta + Zpe * gamma;
        ADREPORT(pred);
        } break;
      default:
        error("Unknown scale!");
    }
  }

  switch (part) {
    case post: {
      nll += npld;
      return nll;
    } break;
    case lik:
      return nll;
      break;
    case pri:
      return npld;
      break;
    default:
      error("Unknown output type!");
  }

  // return nll;
}
