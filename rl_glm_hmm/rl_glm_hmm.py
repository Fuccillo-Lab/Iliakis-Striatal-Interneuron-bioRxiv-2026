#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Hybrid RL-HMM for Iris Stone-style GLM-HMM infrastructure.

Model:
    state-specific alpha controls one Q-learning trace per latent state

    state-specific policy:
        logit_z(t) =
            bias_z
          + beta_pos_z * max(deltaQ_t, 0)
          + beta_neg_z * min(deltaQ_t, 0)
          + sum_j gamma_{z,j} * X_extra[t,j]

This is intended for adding GLM-style history predictors to the RL-HMM,
e.g. prevChoice_prevReward and prevChoice_prevUnreward.

Choice convention:
    y = 1 means push
    y = 0 means pull

Reward convention:
    reward is 0/1 by default; reward_scale=8.0 converts to uL units.

Important options:
    fit_alpha=True:
        theta = [raw_alpha_1..raw_alpha_K, bias, beta_pos, beta_neg, gamma_extra...]

    fit_alpha=False:
        theta = [bias, beta_pos, beta_neg, gamma_extra...]
        alpha = alpha_fixed
"""

import numpy as np
from scipy import optimize

try:
    from glmhmm.hmm import HMM
    from glmhmm.init_params import init_transitions, init_states
except ImportError:
    from hmm import HMM
    from init_params import init_transitions, init_states


def sigmoid(x):
    x = np.clip(x, -50, 50)
    return 1.0 / (1.0 + np.exp(-x))


def logit(p):
    p = np.clip(p, 1e-6, 1 - 1e-6)
    return np.log(p / (1 - p))


# Optional fast path. If numba is unavailable, the class falls back to pure NumPy/Python.
try:
    from numba import njit
    NUMBA_AVAILABLE = True
except Exception:
    NUMBA_AVAILABLE = False
    def njit(*args, **kwargs):
        def deco(fn):
            return fn
        return deco


@njit(cache=True)
def _sigmoid_numba(x):
    if x > 50.0:
        x = 50.0
    elif x < -50.0:
        x = -50.0
    return 1.0 / (1.0 + np.exp(-x))


@njit(cache=True)
def _compute_phi_numba(y, rewards, sess, alpha, bias, beta_pos, beta_neg, gamma_extra,
                       reward_scale, q_init, eps):
    n = y.shape[0]
    k_states = alpha.shape[0]
    n_extra = gamma_extra.shape[1]
    phi = np.zeros((n, k_states, 2), dtype=np.float64)

    for k in range(k_states):
        for si in range(sess.shape[0] - 1):
            s0 = sess[si]
            s1 = sess[si + 1]
            q_pull = q_init
            q_push = q_init

            for t in range(s0, s1):
                delta_q = q_push - q_pull
                dq_pos = delta_q if delta_q > 0.0 else 0.0
                dq_neg = delta_q if delta_q < 0.0 else 0.0

                logit_val = bias[k] + beta_pos[k] * dq_pos + beta_neg[k] * dq_neg
                for j in range(n_extra):
                    logit_val += gamma_extra[k, j] * rewards[0] * 0.0 + gamma_extra[k, j] * 0.0  # overwritten below

                # Separate loop over X was removed from this function signature in an earlier draft.
                # This placeholder should never be used.
                p_push = _sigmoid_numba(logit_val)
                if p_push < eps:
                    p_push = eps
                elif p_push > 1.0 - eps:
                    p_push = 1.0 - eps

                phi[t, k, 1] = p_push
                phi[t, k, 0] = 1.0 - p_push

                r = rewards[t] * reward_scale
                if y[t] == 1:
                    q_push = q_push + alpha[k] * (r - q_push)
                    q_pull = q_pull * (1.0 - alpha[k])
                elif y[t] == 0:
                    q_pull = q_pull + alpha[k] * (r - q_pull)
                    q_push = q_push * (1.0 - alpha[k])
    return phi


@njit(cache=True)
def _compute_phi_numba_X(y, rewards, sess, X_extra, alpha, bias, beta_pos, beta_neg, gamma_extra,
                         reward_scale, q_init, eps):
    n = y.shape[0]
    k_states = alpha.shape[0]
    n_extra = X_extra.shape[1]
    phi = np.zeros((n, k_states, 2), dtype=np.float64)

    for k in range(k_states):
        for si in range(sess.shape[0] - 1):
            s0 = sess[si]
            s1 = sess[si + 1]
            q_pull = q_init
            q_push = q_init

            for t in range(s0, s1):
                delta_q = q_push - q_pull
                dq_pos = delta_q if delta_q > 0.0 else 0.0
                dq_neg = delta_q if delta_q < 0.0 else 0.0

                logit_val = bias[k] + beta_pos[k] * dq_pos + beta_neg[k] * dq_neg
                for j in range(n_extra):
                    logit_val += X_extra[t, j] * gamma_extra[k, j]

                p_push = _sigmoid_numba(logit_val)
                if p_push < eps:
                    p_push = eps
                elif p_push > 1.0 - eps:
                    p_push = 1.0 - eps

                phi[t, k, 1] = p_push
                phi[t, k, 0] = 1.0 - p_push

                r = rewards[t] * reward_scale
                if y[t] == 1:
                    q_push = q_push + alpha[k] * (r - q_push)
                    q_pull = q_pull * (1.0 - alpha[k])
                elif y[t] == 0:
                    q_pull = q_pull + alpha[k] * (r - q_pull)
                    q_push = q_push * (1.0 - alpha[k])
    return phi


@njit(cache=True)
def _neg_expected_loglike_numba(y, rewards, sess, X_extra, gammas, alpha, bias, beta_pos, beta_neg,
                                gamma_extra, reward_scale, q_init, eps):
    n = y.shape[0]
    k_states = alpha.shape[0]
    n_extra = X_extra.shape[1]
    val = 0.0

    for k in range(k_states):
        for si in range(sess.shape[0] - 1):
            s0 = sess[si]
            s1 = sess[si + 1]
            q_pull = q_init
            q_push = q_init

            for t in range(s0, s1):
                delta_q = q_push - q_pull
                dq_pos = delta_q if delta_q > 0.0 else 0.0
                dq_neg = delta_q if delta_q < 0.0 else 0.0

                logit_val = bias[k] + beta_pos[k] * dq_pos + beta_neg[k] * dq_neg
                for j in range(n_extra):
                    logit_val += X_extra[t, j] * gamma_extra[k, j]

                p_push = _sigmoid_numba(logit_val)
                if p_push < eps:
                    p_push = eps
                elif p_push > 1.0 - eps:
                    p_push = 1.0 - eps

                if y[t] == 1:
                    p_obs = p_push
                else:
                    p_obs = 1.0 - p_push
                    if p_obs < eps:
                        p_obs = eps

                val -= gammas[t, k] * np.log(p_obs)

                r = rewards[t] * reward_scale
                if y[t] == 1:
                    q_push = q_push + alpha[k] * (r - q_push)
                    q_pull = q_pull * (1.0 - alpha[k])
                elif y[t] == 0:
                    q_pull = q_pull + alpha[k] * (r - q_pull)
                    q_push = q_push * (1.0 - alpha[k])
    return val


@njit(cache=True)
def _update_transitions_sessionwise_numba(y, alpha_fwd, beta_bwd, cs, A, phi, sess, eps):
    k_states = A.shape[0]
    numer = np.zeros((k_states, k_states), dtype=np.float64)
    denom = np.zeros(k_states, dtype=np.float64)

    for si in range(sess.shape[0] - 1):
        s0 = sess[si]
        s1 = sess[si + 1]
        for t in range(s0, s1 - 1):
            cst = cs[t + 1]
            if cst < eps:
                cst = eps
            for i in range(k_states):
                row_sum = 0.0
                for j in range(k_states):
                    beta_phi = beta_bwd[t + 1, j] * phi[t + 1, j, y[t + 1]]
                    xi = alpha_fwd[t, i] * beta_phi * A[i, j] / cst
                    numer[i, j] += xi
                    row_sum += xi
                denom[i] += row_sum

    A_new = np.zeros((k_states, k_states), dtype=np.float64)
    for i in range(k_states):
        d = denom[i]
        if d < eps:
            d = eps
        rowsum = 0.0
        for j in range(k_states):
            A_new[i, j] = numer[i, j] / d
            rowsum += A_new[i, j]
        if rowsum < eps:
            rowsum = eps
        for j in range(k_states):
            A_new[i, j] /= rowsum
    return A_new


class RLHMMHybrid(HMM):
    """
    Hidden Markov model with RL-derived and GLM-style hybrid emissions.

    Binary classes:
        class 0 = pull
        class 1 = push
    """

    def __init__(
        self,
        n,
        c,
        k,
        n_extra=0,
        reward_scale=8.0,
        q_init=4.0,
        eps=1e-12,
        fit_alpha=True,
        alpha_fixed=0.5,
        split_beta=True,
        verbose=False,
    ):
        super().__init__(n=n, d=0, c=c, k=k)

        if c != 2:
            raise ValueError("RLHMMHybrid currently assumes binary choices with c=2.")

        self.n_extra = int(n_extra)
        self.reward_scale = float(reward_scale)
        self.q_init = float(q_init)
        self.eps = float(eps)
        self.fit_alpha = bool(fit_alpha)
        self.alpha_fixed = float(alpha_fixed)
        self.split_beta = bool(split_beta)
        self.verbose = bool(verbose)

        if self.n_extra < 0:
            raise ValueError("n_extra must be >= 0.")
        if not (0.0 < self.alpha_fixed < 1.0):
            raise ValueError("alpha_fixed must be strictly between 0 and 1.")

    def theta_length(self):
        n = 0
        if self.fit_alpha:
            n += self.k                    # state-specific raw alphas
        n += self.k                          # bias
        n += (2 * self.k if self.split_beta else self.k)  # beta_pos/beta_neg OR single beta
        n += self.k * self.n_extra           # state-specific extra GLM weights
        return n

    def generate_params(self, transitions=['dirichlet', 10, 1], state_priors='uniform', seed_scale=0.2):
        A = init_transitions(
            self,
            distribution=transitions[0],
            alpha_diag=transitions[1],
            alpha_full=transitions[2],
        )
        pi0 = init_states(self, state_priors)

        theta = np.zeros(self.theta_length(), dtype=float)
        ix = 0

        if self.fit_alpha:
            theta[ix:ix + self.k] = logit(self.alpha_fixed)
            ix += self.k

        # bias, beta(s), and extras
        theta[ix:] = np.random.uniform(-seed_scale, seed_scale, size=theta.size - ix)

        return A, theta, pi0

    def unpack_theta(self, theta):
        """
        Return:
            alpha
            bias             shape (K,)
            beta_pos         shape (K,)
            beta_neg         shape (K,)
            gamma_extra      shape (K, n_extra)
        """
        theta = np.asarray(theta, dtype=float)

        expected = self.theta_length()
        if theta.size != expected:
            raise ValueError(f"Expected theta length {expected}, got {theta.size}.")

        ix = 0

        if self.fit_alpha:
            alpha = sigmoid(theta[ix:ix + self.k])
            ix += self.k
        else:
            alpha = np.full(self.k, self.alpha_fixed, dtype=float)

        bias = theta[ix:ix + self.k]
        ix += self.k

        if self.split_beta:
            beta_pos = theta[ix:ix + self.k]
            ix += self.k
            beta_neg = theta[ix:ix + self.k]
            ix += self.k
        else:
            beta_single = theta[ix:ix + self.k]
            ix += self.k
            beta_pos = beta_single
            beta_neg = beta_single

        if self.n_extra > 0:
            gamma_extra = theta[ix:ix + self.k * self.n_extra].reshape(self.k, self.n_extra)
        else:
            gamma_extra = np.zeros((self.k, 0), dtype=float)

        return alpha, bias, beta_pos, beta_neg, gamma_extra

    def compute_q_trace(self, y, rewards, sess, alpha):
        y = np.asarray(y).astype(int)
        rewards = np.asarray(rewards, dtype=float) * self.reward_scale
        Q = np.zeros((len(y), 2), dtype=float)

        for s0, s1 in zip(sess[:-1], sess[1:]):
            q_pull = self.q_init
            q_push = self.q_init

            for t in range(s0, s1):
                Q[t, 0] = q_pull
                Q[t, 1] = q_push

                r = rewards[t]

                if y[t] == 1:      # push chosen
                    q_push = q_push + alpha * (r - q_push)
                    q_pull = q_pull * (1.0 - alpha)
                elif y[t] == 0:    # pull chosen
                    q_pull = q_pull + alpha * (r - q_pull)
                    q_push = q_push * (1.0 - alpha)
                else:
                    pass

        return Q

    def compute_phi(self, y, rewards, sess, theta, X_extra=None):
        """
        Return phi[t,k,c] = P(choice class c at trial t | latent state k).
        Fast path uses numba when installed.
        """
        alpha, bias, beta_pos, beta_neg, gamma_extra = self.unpack_theta(theta)
        y = np.asarray(y, dtype=np.int64)
        rewards = np.asarray(rewards, dtype=np.float64)
        sess = np.asarray(sess, dtype=np.int64)

        if self.n_extra > 0:
            if X_extra is None:
                raise ValueError("X_extra required when n_extra > 0.")
            X_extra = np.asarray(X_extra, dtype=np.float64)
            if X_extra.shape != (len(y), self.n_extra):
                raise ValueError(
                    f"Expected X_extra shape {(len(y), self.n_extra)}, got {X_extra.shape}."
                )
        else:
            X_extra = np.zeros((len(y), 0), dtype=np.float64)

        if NUMBA_AVAILABLE:
            return _compute_phi_numba_X(
                y, rewards, sess, X_extra,
                np.asarray(alpha, dtype=np.float64),
                np.asarray(bias, dtype=np.float64),
                np.asarray(beta_pos, dtype=np.float64),
                np.asarray(beta_neg, dtype=np.float64),
                np.asarray(gamma_extra, dtype=np.float64),
                self.reward_scale, self.q_init, self.eps,
            )

        # Pure Python fallback.
        phi = np.zeros((len(y), self.k, self.c), dtype=float)
        for k in range(self.k):
            Qk = self.compute_q_trace(y, rewards, sess, alpha[k])
            delta_q = Qk[:, 1] - Qk[:, 0]
            dq_pos = np.maximum(delta_q, 0.0)
            dq_neg = np.minimum(delta_q, 0.0)
            logits = bias[k] + beta_pos[k] * dq_pos + beta_neg[k] * dq_neg
            if self.n_extra > 0:
                logits = logits + X_extra @ gamma_extra[k, :]
            p_push = sigmoid(logits)
            p_push = np.clip(p_push, self.eps, 1.0 - self.eps)
            phi[:, k, 1] = p_push
            phi[:, k, 0] = 1.0 - p_push
        return phi

    def neg_expected_loglike(self, theta, y, rewards, sess, gammas, X_extra=None, l2=0.0, l2_extra=0.0):
        alpha, bias, beta_pos, beta_neg, gamma_extra = self.unpack_theta(theta)
        y = np.asarray(y, dtype=np.int64)
        rewards = np.asarray(rewards, dtype=np.float64)
        sess = np.asarray(sess, dtype=np.int64)
        gammas = np.asarray(gammas, dtype=np.float64)

        if self.n_extra > 0:
            X_extra = np.asarray(X_extra, dtype=np.float64)
        else:
            X_extra = np.zeros((len(y), 0), dtype=np.float64)

        if NUMBA_AVAILABLE:
            val = _neg_expected_loglike_numba(
                y, rewards, sess, X_extra, gammas,
                np.asarray(alpha, dtype=np.float64),
                np.asarray(bias, dtype=np.float64),
                np.asarray(beta_pos, dtype=np.float64),
                np.asarray(beta_neg, dtype=np.float64),
                np.asarray(gamma_extra, dtype=np.float64),
                self.reward_scale, self.q_init, self.eps,
            )
        else:
            phi = self.compute_phi(y, rewards, sess, theta, X_extra=X_extra)
            idx = y.astype(int)
            p_obs = phi[np.arange(len(y)), :, idx]
            val = -np.sum(gammas * np.log(np.clip(p_obs, self.eps, 1.0)))

        if l2 > 0:
            start = self.k if self.fit_alpha else 0
            n_base = self.k + (2 * self.k if self.split_beta else self.k)
            end = start + n_base
            val = val + 0.5 * l2 * np.sum(theta[start:end] ** 2)

        if self.n_extra > 0 and l2_extra > 0:
            n_base = self.k + (2 * self.k if self.split_beta else self.k)
            start = (self.k if self.fit_alpha else 0) + n_base
            val = val + 0.5 * l2_extra * np.sum(theta[start:] ** 2)

        return float(val)

    def update_observations(self, y, rewards, sess, theta, gammas, X_extra=None,
                            l2=0.0, l2_extra=0.0, maxiter=500):
        bounds = []

        if self.fit_alpha:
            # Bound alpha itself to [0.05, 0.95].
            # alpha = sigmoid(raw_alpha), so convert alpha bounds to raw-logit bounds.
            alpha_min = 0.05
            alpha_max = 0.95
            bounds += [(logit(alpha_min), logit(alpha_max))] * self.k

        n_base = self.k + (2 * self.k if self.split_beta else self.k)
        bounds += [(-20, 20)] * n_base

        if self.n_extra > 0:
            bounds += [(-20, 20)] * (self.k * self.n_extra)

        if self.verbose:
            alpha, _, _, _, gamma_extra = self.unpack_theta(theta)
            print(
                f"    obs M-step start: fit_alpha={self.fit_alpha}, "
                f"alpha={np.round(alpha, 6)}, n_extra={self.n_extra}"
            )

        result = optimize.minimize(
            fun=lambda th: self.neg_expected_loglike(
                th, y, rewards, sess, gammas,
                X_extra=X_extra,
                l2=l2,
                l2_extra=l2_extra,
            ),
            x0=np.asarray(theta, dtype=float),
            method='L-BFGS-B',
            bounds=bounds,
            options={'maxiter': maxiter, 'ftol': 1e-9, 'maxls': 50},
        )

        theta_new = result.x
        phi_new = self.compute_phi(y, rewards, sess, theta_new, X_extra=X_extra)
        self.last_obs_opt_result = result

        if self.verbose:
            alpha, _, _, _, _ = self.unpack_theta(theta_new)
            print(
                f"    obs M-step done : success={result.success}, "
                f"fun={result.fun:.6f}, alpha={np.round(alpha, 6)}"
            )

        return theta_new, phi_new

    def update_transitions_sessionwise(self, y, alpha_fwd, beta_bwd, cs, A, phi, sess):
        if NUMBA_AVAILABLE:
            return _update_transitions_sessionwise_numba(
                np.asarray(y, dtype=np.int64),
                np.asarray(alpha_fwd, dtype=np.float64),
                np.asarray(beta_bwd, dtype=np.float64),
                np.asarray(cs, dtype=np.float64),
                np.asarray(A, dtype=np.float64),
                np.asarray(phi, dtype=np.float64),
                np.asarray(sess, dtype=np.int64),
                self.eps,
            )

        numer = np.zeros((self.k, self.k), dtype=float)
        denom = np.zeros((self.k, 1), dtype=float)

        for s0, s1 in zip(sess[:-1], sess[1:]):
            for t in range(s0, s1 - 1):
                beta_phi = beta_bwd[t + 1, :] * phi[t + 1, :, int(y[t + 1])]
                xi = (
                    (alpha_fwd[t, :].reshape(self.k, 1) * beta_phi.reshape(1, self.k))
                    * A
                ) / np.clip(cs[t + 1], self.eps, np.inf)

                numer += xi
                denom += np.sum(xi, axis=1, keepdims=True)

        A_new = numer / np.clip(denom, self.eps, np.inf)
        A_new = A_new / np.sum(A_new, axis=1, keepdims=True)
        return A_new

    def fit(self, y, rewards, A, theta, pi0=None, fit_init_states=False,
            maxiter=100, tol=1e-3, sess=None, B=1, X_extra=None,
            l2=0.0, l2_extra=0.0, obs_maxiter=500):
        y = np.asarray(y).astype(int)
        rewards = np.asarray(rewards, dtype=float)

        if sess is None:
            sess = np.array([0, len(y)], dtype=int)
        else:
            sess = np.asarray(sess, dtype=int)

        if self.n_extra > 0:
            X_extra = np.asarray(X_extra, dtype=float)
            if X_extra.shape != (len(y), self.n_extra):
                raise ValueError(
                    f"Expected X_extra shape {(len(y), self.n_extra)}, got {X_extra.shape}."
                )
        else:
            X_extra = None

        self.lls = np.empty(maxiter)
        self.lls[:] = np.nan
        self.pi0 = pi0

        if self.verbose:
            print(
                f"RLHMMHybrid.fit: fit_alpha={self.fit_alpha}, "
                f"alpha_fixed={self.alpha_fixed}, split_beta={self.split_beta}, n_extra={self.n_extra}"
            )

        phi = self.compute_phi(y, rewards, sess, theta, X_extra=X_extra)

        for it in range(maxiter):
            alpha_fwd = np.zeros((self.n, self.k))
            beta_bwd = np.zeros_like(alpha_fwd)
            cs = np.zeros(self.n)
            self.pStates = np.zeros_like(alpha_fwd)
            self.states = np.zeros(self.n)
            ll = 0.0

            for s in range(len(sess) - 1):
                sl = slice(sess[s], sess[s + 1])
                ll_s, alpha_s, _, cs_s = self.forwardPass(y[sl], A, phi[sl, :, :], pi0=pi0)
                pback_s, beta_s, zhat_s = self.backwardPass(
                    y[sl], A, phi[sl, :, :], alpha_s, cs_s
                )

                ll += ll_s
                alpha_fwd[sl] = alpha_s
                beta_bwd[sl] = beta_s
                cs[sl] = cs_s
                self.pStates[sl] = pback_s ** B
                self.states[sl] = zhat_s

            self.lls[it] = ll

            if self.verbose:
                alpha_curr, _, _, _, _ = self.unpack_theta(theta)
                occ = self.pStates.mean(axis=0)
                print(
    f"  EM {it:03d}: LL={ll:.6f}, "
    f"alpha={np.round(alpha_curr, 4)}, "
    f"occ={np.round(occ, 4)}"
)

            A = self.update_transitions_sessionwise(y, alpha_fwd, beta_bwd, cs, A, phi, sess)
            theta, phi = self.update_observations(
                y,
                rewards,
                sess,
                theta,
                self.pStates,
                X_extra=X_extra,
                l2=l2,
                l2_extra=l2_extra,
                maxiter=obs_maxiter,
            )

            if fit_init_states:
                pi0 = self._updateInitStates(self.pStates)

            if it > 5 and self.lls[it - 5] + tol >= ll:
                break

        self.A = A
        self.theta = theta
        self.phi = phi
        self.pi0 = pi0
        return self.lls, self.A, self.theta, self.pi0
