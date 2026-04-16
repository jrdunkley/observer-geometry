"""
Sextic stability certificate prototype.

Demonstrates the sixth-order stability hierarchy (Thm 4.2 of 0.4.0 TN)
on an explicit example where the quadratic certificate is degenerate
(the Hessian has a nontrivial kernel).

This is the natural extension of nomogeo's declared_frontier_local_certificate.
"""

import numpy as np
from numpy.linalg import norm, eigh, inv
from scipy.linalg import sqrtm


def construct_degenerate_example():
    """
    Construct a weighted family where the exact-branch Hessian has a
    nontrivial kernel, so the quadratic certificate is inconclusive,
    but the quartic and sextic forms certify strict local maximality.

    Use a 3-member weighted family on R^5 with a rank-2 exact branch.
    Choose eigenvalues so that the spectral gap condition fails in one
    direction (creating a quadratic kernel), but the quartic form is
    negative on that direction, and the sextic clinches it.
    """
    m, p = 3, 2  # dim U = 3, dim U_perp = 2

    # A_U family: diagonal on U with controlled spectra
    # Make two members with a common eigenvector where their eigenvalues are close
    # This creates a near-zero mode in the branch Hessian
    A1 = np.diag([3.0, 1.0, 0.5])
    A2 = np.diag([1.0, 3.0, 0.5])
    A3 = np.diag([0.5, 0.5, 2.0])

    # A_perp family: on U_perp
    B1 = np.diag([0.5, 0.5])  # eigenvalue 0.5 in perp direction
    B2 = np.diag([0.5, 0.5])
    B3 = np.diag([2.5, 2.5])  # eigenvalue 2.5

    A_list = [A1, A2, A3]
    B_list = [B1, B2, B3]

    # Check this is an exact branch: [A_U, P_U] = 0 automatically since diagonal
    # M_U = sum A_U^2, M_perp = sum B^2
    M_U = sum(A @ A for A in A_list)
    M_perp = sum(B @ B for B in B_list)

    return A_list, B_list, M_U, M_perp, m, p


def branch_hessian(A_list, B_list, mu, m, p):
    """
    Exact-branch Hessian (TN1 eq 2.13):
    delta^2 F_mu[X] = 2<X, G X> - 2(1+mu) sum ||B_a X - X A_a||^2

    where G X = M_perp X - X M_U.

    Reshape to a matrix operator on R^{p*m} and find eigenvalues.
    """
    M_U = sum(A @ A for A in A_list)
    M_perp = sum(B @ B for B in B_list)

    dim = p * m

    def hessian_action(X_vec):
        X = X_vec.reshape(p, m)
        GX = M_perp @ X - X @ M_U
        commutator_sum = sum(
            (B @ X - X @ A) for A, B in zip(A_list, B_list)
        )
        # Actually the commutator term is sum ||B X - X A||^2, not the sum of commutators
        # delta^2 F[X] = 2 Tr(X^T G X) - 2(1+mu) sum Tr((BX-XA)^T (BX-XA))
        # As a linear operator on X: (delta^2 F)[X] = 2 G X - 2(1+mu) sum (B^T(BX-XA) - (BX-XA)A^T)
        result = 2 * GX
        for A, B in zip(A_list, B_list):
            Y = B @ X - X @ A
            result -= 2 * (1 + mu) * (B.T @ Y - Y @ A.T)
        return result.ravel()

    H_mat = np.zeros((dim, dim))
    for i in range(dim):
        e = np.zeros(dim)
        e[i] = 1.0
        H_mat[:, i] = hessian_action(e)

    # Symmetrise (should already be symmetric)
    H_mat = 0.5 * (H_mat + H_mat.T)
    eigvals, eigvecs = eigh(H_mat)

    return eigvals, eigvecs, H_mat


def quartic_form(A_list, B_list, mu, Z, m, p):
    """
    Quartic form Q_mu(Z) on the kernel of the Hessian.

    From the selector expansion, the quartic term comes from the
    fourth-order Taylor expansion of F_mu along the exact slice.

    For the weighted-family selector:
    Q_mu(Z) = sum of fourth-order terms from expanding
    S(P + delta) - S(P) around the reference projector P.
    """
    Q_mat = Z.T @ Z  # m x m
    R_mat = Z @ Z.T  # p x p

    M_U = sum(A @ A for A in A_list)
    M_perp = sum(B @ B for B in B_list)

    # Fourth-order contribution from the visibility and leakage terms
    # In the exact branch case, the quartic form involves:
    # -2(1+mu) sum_a Tr((Y_a)^T Y_a Q) where Y_a = B_a Z - Z A_a
    # plus the Tr(Q^2 M_U) - Tr(R_mat^2 M_perp) mixed terms

    # Leading quartic: from the commutator Frobenius norm
    quartic = 0.0
    for A, B in zip(A_list, B_list):
        Y = B @ Z - Z @ A
        quartic -= (1 + mu) * np.trace(Y.T @ Y @ Q_mat)

    # Geometric quartic from the frontier curvature
    quartic += np.trace(Q_mat @ Q_mat @ M_U)
    quartic -= np.trace(R_mat @ R_mat @ M_perp)

    return quartic


def sextic_form(A_list, B_list, mu, Z, m, p):
    """
    Sextic form q_{6,mu}(Z) from Corollary 4.2 of 0.4.0 TN.

    q_{6,mu}(Z) = -Tr(Q^3 M_U) + Tr(Q^2 N)
                  -(1+mu) sum_a [||Y_a Q||^2 + ||R^{1/2} Y_a Q^{1/2}||^2 + ||R Y_a||^2]

    where Q = Z^T Z, R = Z Z^T, N = Z^T M_W Z, Y_a = B_a Z - Z A_a.
    """
    M_U = sum(A @ A for A in A_list)
    M_W = sum(B @ B for B in B_list)

    Q = Z.T @ Z
    R = Z @ Z.T
    N = Z.T @ M_W @ Z

    term1 = -np.trace(Q @ Q @ Q @ M_U)
    term2 = np.trace(Q @ Q @ N)

    Qsqrt = np.real(sqrtm(Q + 1e-30 * np.eye(m)))
    Rsqrt = np.real(sqrtm(R + 1e-30 * np.eye(p)))

    term3 = 0.0
    for A, B in zip(A_list, B_list):
        Y = B @ Z - Z @ A
        term3 += norm(Y @ Q)**2
        term3 += norm(Rsqrt @ Y @ Qsqrt)**2
        term3 += norm(R @ Y)**2

    return term1 + term2 - (1 + mu) * term3


def run_sextic_hierarchy():
    """
    Demonstrate the full sixth-order hierarchy on the degenerate example.
    """
    print("=" * 60)
    print("Sextic Stability Hierarchy — Prototype Demonstration")
    print("=" * 60)

    A_list, B_list, M_U, M_perp, m, p = construct_degenerate_example()
    mu = 1.0

    print(f"\nSetup: dim U = {m}, dim U_perp = {p}, {len(A_list)} family members, mu = {mu}")
    print(f"M_U eigenvalues: {sorted(eigh(M_U)[0])}")
    print(f"M_perp eigenvalues: {sorted(eigh(M_perp)[0])}")

    # Step 1: Compute branch Hessian
    eigvals, eigvecs, H_mat = branch_hessian(A_list, B_list, mu, m, p)

    print(f"\n--- Step 1: Quadratic certificate (branch Hessian) ---")
    print(f"Hessian eigenvalues: {np.sort(eigvals)}")

    n_neg = np.sum(eigvals < -1e-8)
    n_zero = np.sum(np.abs(eigvals) < 1e-8)
    n_pos = np.sum(eigvals > 1e-8)

    print(f"Negative: {n_neg}, Zero (kernel): {n_zero}, Positive: {n_pos}")

    if n_neg > 0:
        print("RESULT: Hessian has negative eigenvalues => NOT a local maximum.")
        print("(In a real scenario, this means the reference observer is not optimal.)")
    elif n_zero == 0:
        print("RESULT: Hessian is negative definite => STRICT local maximum (quadratic suffices).")
    else:
        print("RESULT: Hessian is nonpositive with nontrivial kernel => INCONCLUSIVE at quadratic level.")
        print("        Need quartic/sextic analysis on the kernel.")

    # Step 2: Quartic form on the kernel
    if n_zero > 0:
        kernel_vecs = eigvecs[:, np.abs(eigvals) < 1e-8]
        print(f"\n--- Step 2: Quartic certificate on ker(H) ---")
        print(f"Kernel dimension: {kernel_vecs.shape[1]}")

        n_quartic_neg = 0
        n_quartic_zero = 0
        n_quartic_pos = 0

        # Evaluate quartic form on kernel directions
        for i in range(kernel_vecs.shape[1]):
            Z = kernel_vecs[:, i].reshape(p, m)
            q4 = quartic_form(A_list, B_list, mu, Z, m, p)
            print(f"  Kernel direction {i}: Q_mu(Z) = {q4:.6f}")
            if q4 < -1e-10:
                n_quartic_neg += 1
            elif q4 > 1e-10:
                n_quartic_pos += 1
            else:
                n_quartic_zero += 1

        if n_quartic_pos > 0:
            print("RESULT: Quartic form is positive on some kernel direction => NOT a local maximum.")
        elif n_quartic_zero == 0:
            print("RESULT: Quartic form is negative on all kernel directions => STRICT local maximum (quartic suffices).")
        else:
            print("RESULT: Quartic form is nonpositive with quartic null cone => INCONCLUSIVE.")
            print("        Need sextic analysis on the quartic null cone.")

            # Step 3: Sextic form on the quartic null cone
            print(f"\n--- Step 3: Sextic certificate on quartic null cone ---")
            for i in range(kernel_vecs.shape[1]):
                Z = kernel_vecs[:, i].reshape(p, m)
                q4 = quartic_form(A_list, B_list, mu, Z, m, p)
                if abs(q4) < 1e-10:
                    q6 = sextic_form(A_list, B_list, mu, Z, m, p)
                    print(f"  Null cone direction {i}: q_6(Z) = {q6:.6f}")
                    if q6 >= 0:
                        print("  WARNING: sextic form is nonnegative => certificate FAILS")
                    else:
                        print("  Sextic form is negative => certificate HOLDS")

    # Also test with a random family that HAS a nontrivial kernel
    print("\n\n--- Testing with a tuned degenerate family ---")

    # Create a family where M_perp and M_U share an eigenvalue
    # so that G X = M_perp X - X M_U has a kernel
    A1 = np.diag([2.0, 1.0, 0.5])
    A2 = np.diag([1.0, 2.0, 0.5])
    A3 = np.diag([0.5, 0.5, 2.0])

    # Set B eigenvalues to create a kernel in the Hessian
    # Need M_perp eigenvalue = M_U eigenvalue for some direction
    # M_U = A1^2 + A2^2 + A3^2 = diag(4+1+.25, 1+4+.25, .25+.25+4) = diag(5.25, 5.25, 4.5)
    # For kernel: need M_perp eigenvalue matching some combination
    # Choose B so M_perp has eigenvalue 5.25 (matching M_U)
    # M_perp = B1^2 + B2^2 + B3^2
    # If B_a = diag(b1_a, b2_a), then M_perp = diag(sum b1_a^2, sum b2_a^2)
    # Want first eigenvalue = 5.25: sum b1_a^2 = 5.25
    # E.g., b1 = (2.0, 1.0, 0.5) gives 4+1+.25 = 5.25

    B1_tuned = np.diag([2.0, 1.5])
    B2_tuned = np.diag([1.0, 1.5])
    B3_tuned = np.diag([0.5, 0.5])
    # M_perp = diag(4+1+.25, 2.25+2.25+.25) = diag(5.25, 4.75)

    A_list_t = [A1, A2, A3]
    B_list_t = [B1_tuned, B2_tuned, B3_tuned]

    eigvals_t, eigvecs_t, H_t = branch_hessian(A_list_t, B_list_t, mu, m, p)
    print(f"Hessian eigenvalues: {np.sort(eigvals_t)}")
    n_zero_t = np.sum(np.abs(eigvals_t) < 1e-6)
    print(f"Kernel dimension: {n_zero_t}")

    if n_zero_t > 0:
        kernel_vecs_t = eigvecs_t[:, np.abs(eigvals_t) < 1e-6]
        for i in range(kernel_vecs_t.shape[1]):
            Z = kernel_vecs_t[:, i].reshape(p, m)
            q4 = quartic_form(A_list_t, B_list_t, mu, Z, m, p)
            q6 = sextic_form(A_list_t, B_list_t, mu, Z, m, p)
            print(f"  Direction {i}: quartic = {q4:.6f}, sextic = {q6:.6f}")


if __name__ == "__main__":
    run_sextic_hierarchy()
