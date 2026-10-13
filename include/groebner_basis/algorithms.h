#pragma once

#include <stddef.h>

#include <algorithm>
#include <cassert>
#include <ostream>
#include <set>
#include <utility>
#include <vector>

#include "polynomial.h"
#include "term.h"

namespace groebner_basis::algorithms {

template <typename Coefficient, typename MonomialOrder>
using PolynomialBasis = std::vector<Polynomial<Coefficient, MonomialOrder>>;

template <typename Coefficient, typename MonomialOrder>
Polynomial<Coefficient, MonomialOrder> ReducePolynomial(
    Polynomial<Coefficient, MonomialOrder> dividend,
    const std::vector<Polynomial<Coefficient, MonomialOrder>>& divisors
);

template <typename Coefficient, typename MonomialOrder>
Polynomial<Coefficient, MonomialOrder> ComputeSPolynomial(
    const Polynomial<Coefficient, MonomialOrder>& lhs,
    const Polynomial<Coefficient, MonomialOrder>& rhs
);

namespace detail {

using BuchbergerPair = std::pair<size_t, size_t>;
using BuchbergerPairSet = std::set<BuchbergerPair>;

inline BuchbergerPair MakeBuchbergerPair(size_t lhs, size_t rhs) {
    assert(lhs != rhs);

    if (lhs < rhs) {
        return {lhs, rhs};
    }
    return {rhs, lhs};
}

template <typename Coefficient, typename MonomialOrder>
BuchbergerPairSet BuildBuchbergerPairs(const PolynomialBasis<Coefficient, MonomialOrder>& basis) {
    BuchbergerPairSet result;

    for (size_t lhs = 0; lhs != basis.size(); ++lhs) {
        for (size_t rhs = lhs + 1; rhs != basis.size(); ++rhs) {
            result.emplace_hint(result.end(), lhs, rhs);
        }
    }

    return result; 
}

template <typename Coefficient, typename MonomialOrder>
bool MeetsBuchbergerProductCriterion(
    const Polynomial<Coefficient, MonomialOrder>& lhs,
    const Polynomial<Coefficient, MonomialOrder>& rhs
) {
    const auto& lhs_leading_monomial = lhs.GetLeadingMonomial();
    const auto& rhs_leading_monomial = rhs.GetLeadingMonomial();

    const auto least_common_multiple =
        ComputeLeastCommonMultiple(lhs_leading_monomial, rhs_leading_monomial);

    return least_common_multiple == lhs_leading_monomial * rhs_leading_monomial;
}

template <typename Coefficient, typename MonomialOrder>
bool MeetsBuchbergerChainCriterion(
    size_t lhs_idx,
    size_t rhs_idx,
    const PolynomialBasis<Coefficient, MonomialOrder>& basis,
    const BuchbergerPairSet& pending_pairs
) {
    const auto least_common_multiple =
        ComputeLeastCommonMultiple(
            basis[lhs_idx].GetLeadingMonomial(),
            basis[rhs_idx].GetLeadingMonomial()
        );

    for (size_t middle_idx = 0; middle_idx != basis.size(); ++middle_idx) {
        if (middle_idx == lhs_idx || middle_idx == rhs_idx) {
            continue;
        }

        const auto lhs_middle_pair = MakeBuchbergerPair(lhs_idx, middle_idx);
        const auto rhs_middle_pair = MakeBuchbergerPair(rhs_idx, middle_idx);

        if (pending_pairs.contains(lhs_middle_pair) || pending_pairs.contains(rhs_middle_pair)) {
            continue;
        }

        if (basis[middle_idx].GetLeadingMonomial().Divides(least_common_multiple)) {
            return true;
        }
    }

    return false;
}

enum class Improvement : unsigned char {
    kNone,
    kProductCriterion,
    kProductAndChainCriteria,
};

template <typename Coefficient, typename MonomialOrder>
PolynomialBasis<Coefficient, MonomialOrder> RunBuchbergerAlgorithm(
    PolynomialBasis<Coefficient, MonomialOrder> generators,
    Improvement improvement = Improvement::kProductAndChainCriteria
) {
    const size_t generator_count = generators.size();
    auto basis = std::move(generators);
    basis.erase(
        std::remove_if(
            basis.begin(), basis.end(),
            [](const auto& polynomial) {
                return polynomial.IsZero();
            }
        ),
        basis.end()
    );

    auto pending_pairs = BuildBuchbergerPairs(basis);

    while (!pending_pairs.empty()) {
        const auto pair_position = pending_pairs.begin();
        const auto [lhs_idx, rhs_idx] = *pair_position;
        pending_pairs.erase(pair_position);

        if (improvement != Improvement::kNone &&
            MeetsBuchbergerProductCriterion(basis[lhs_idx], basis[rhs_idx])) {
            continue;
        }

        if (improvement == Improvement::kProductAndChainCriteria &&
            MeetsBuchbergerChainCriterion(lhs_idx, rhs_idx, basis, pending_pairs)) {
            continue;
        }

        auto remainder = ReducePolynomial(
            ComputeSPolynomial(basis[lhs_idx], basis[rhs_idx]), basis
        );

        if (remainder.IsZero()) {
            continue;
        }

        const size_t new_idx = basis.size();
        basis.push_back(std::move(remainder));
        for (size_t existing_idx = 0; existing_idx != new_idx; ++existing_idx) {
            pending_pairs.emplace(existing_idx, new_idx);
        }
    }

    return basis;
}

template <typename Coefficient, typename MonomialOrder>
Polynomial<Coefficient, MonomialOrder> ReducePolynomialImpl(
    Polynomial<Coefficient, MonomialOrder> dividend,
    const PolynomialBasis<Coefficient, MonomialOrder>& divisors,
    const Polynomial<Coefficient, MonomialOrder>* excluded_divisor = nullptr
) {
    Polynomial<Coefficient, MonomialOrder> remainder;

    while (!dividend.IsZero()) {
        const auto leading_term = dividend.GetLeadingTerm();

        bool reduced = false;
        for (const auto& divisor : divisors) {
            if (&divisor == excluded_divisor || divisor.IsZero()) {
                continue;
            }

            const auto divisor_leading_term = divisor.GetLeadingTerm();
            if (divisor_leading_term.Divides(leading_term)) {
                dividend -= (leading_term / divisor_leading_term) * divisor;
                reduced = true;
                break;
            }
        }

        if (!reduced) {
            remainder += leading_term;
            dividend -= leading_term;
        }
    }

    return remainder;
}

template <typename Coefficient, typename MonomialOrder>
PolynomialBasis<Coefficient, MonomialOrder> MinimizeGroebnerBasis(
    PolynomialBasis<Coefficient, MonomialOrder> basis
) {
    std::vector<bool> keep(basis.size(), true);
    for (size_t candidate_idx = 0; candidate_idx != basis.size(); ++candidate_idx) {
        const auto& candidate_lm = basis[candidate_idx].GetLeadingMonomial();

        bool is_redundant = false;
        for (size_t current_idx = 0; current_idx != basis.size(); ++current_idx) {
            if (candidate_idx == current_idx) {
                continue;
            }

            const auto& current_lm = basis[current_idx].GetLeadingMonomial();

            const bool divides_candidate = current_lm.Divides(candidate_lm);
            const bool takes_precedence = current_lm != candidate_lm || current_idx < candidate_idx;
            if (divides_candidate && takes_precedence) {
                keep[candidate_idx] = false;
                break;
            }
        }
    }

    size_t write_idx = 0;
    for (size_t read_idx = 0; read_idx != basis.size(); ++read_idx) {
        if (!keep[read_idx]) {
            continue;
        }

        if (write_idx != read_idx) {
            basis[write_idx] = std::move(basis[read_idx]);
        }

        ++write_idx;
    }
    basis.resize(write_idx);

    for (auto& polynomial : basis) {
        polynomial.MakeMonic();
        
        if (polynomial.IsOne()) {
            return {std::move(polynomial)};
        }
    } 

    return basis;
}

template <typename Coefficient, typename MonomialOrder>
void SortBasisByLeadingMonomial(PolynomialBasis<Coefficient, MonomialOrder>& basis) {
    std::stable_sort(
        basis.begin(), basis.end(),
        [](const auto& lhs, const auto& rhs) {
            return MonomialOrder{}(
                lhs.GetLeadingMonomial(), rhs.GetLeadingMonomial()
            );
        }
    );
}

template <typename Coefficient, typename MonomialOrder>
PolynomialBasis<Coefficient, MonomialOrder> InterreduceGroebnerBasis(
    PolynomialBasis<Coefficient, MonomialOrder> minimal_basis
) {
    PolynomialBasis<Coefficient, MonomialOrder> result;
    result.reserve(minimal_basis.size());

    for (size_t polynomial_idx = 0; polynomial_idx != minimal_basis.size(); ++polynomial_idx) {
        auto reduced = ReducePolynomialImpl(
            minimal_basis[polynomial_idx], minimal_basis, &minimal_basis[polynomial_idx]
        );

        assert(!reduced.IsZero());
        assert(reduced.IsMonic());

        result.push_back(std::move(reduced));
    }

    SortBasisByLeadingMonomial(result);
    return result;
}

}  // namespace detail

template <typename Coefficient, typename MonomialOrder>
Polynomial<Coefficient, MonomialOrder> ReducePolynomial(
    Polynomial<Coefficient, MonomialOrder> dividend,
    const std::vector<Polynomial<Coefficient, MonomialOrder>>& divisors
) {
    return detail::ReducePolynomialImpl(std::move(dividend), divisors);
}

template <typename Coefficient, typename MonomialOrder>
Polynomial<Coefficient, MonomialOrder> ComputeSPolynomial(
    const Polynomial<Coefficient, MonomialOrder>& lhs,
    const Polynomial<Coefficient, MonomialOrder>& rhs
) {
    assert(!lhs.IsZero() && !rhs.IsZero() && "S-polynomial requires two nonzero polynomials");

    const auto least_common_multiple = ComputeLeastCommonMultiple(
        lhs.GetLeadingMonomial(), rhs.GetLeadingMonomial()
    );
    const Term<Coefficient> common_term{least_common_multiple};

    const auto lhs_factor = common_term / lhs.GetLeadingTerm();
    const auto rhs_factor = common_term / rhs.GetLeadingTerm();
    
    return lhs * lhs_factor - rhs * rhs_factor;
}

template <typename Coefficient, typename MonomialOrder>
PolynomialBasis<Coefficient, MonomialOrder> ComputeGroebnerBasis(
    PolynomialBasis<Coefficient, MonomialOrder> generators
) {
    return detail::RunBuchbergerAlgorithm(std::move(generators));
}

template <typename Coefficient, typename MonomialOrder>
PolynomialBasis<Coefficient, MonomialOrder> ComputeMinimalGroebnerBasis(
    PolynomialBasis<Coefficient, MonomialOrder> generators
) {
    return detail::MinimizeGroebnerBasis(ComputeGroebnerBasis(std::move(generators)));
}

template <typename Coefficient, typename MonomialOrder>
PolynomialBasis<Coefficient, MonomialOrder> ComputeReducedGroebnerBasis(
    PolynomialBasis<Coefficient, MonomialOrder> generators
) {
    return detail::InterreduceGroebnerBasis(ComputeMinimalGroebnerBasis(std::move(generators)));
}

}  // namespace groebner_basis::algorithms
