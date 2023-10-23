/**
 * tlibs2
 * magnon dynamics
 * @author Tobias Weber <tweber@ill.fr>
 * @date january-2022
 * @license GPLv3, see 'LICENSE' file
 *
 * References:
 *   - (Toth 2015) S. Toth and B. Lake, J. Phys.: Condens. Matter 27 166002 (2015):
 *                 https://doi.org/10.1088/0953-8984/27/16/166002
 *                 https://arxiv.org/abs/1402.6069
 *   - (Heinsdorf 2021) N. Heinsdorf, manual example calculation for a simple
 *                      ferromagnetic case, personal communications, 2021/2022.
 *
 * @desc This file implements the formalism given by (Toth 2015).
 * @desc For further references, see the 'LITERATURE' file.
 *
 * ----------------------------------------------------------------------------
 * tlibs
 * Copyright (C) 2017-2023  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * Copyright (C) 2015-2017  Tobias WEBER (Technische Universitaet Muenchen
 *                          (TUM), Garching, Germany).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * ----------------------------------------------------------------------------
 */

#ifndef __TLIBS2_MAGDYN_H__
#define __TLIBS2_MAGDYN_H__

#include <vector>
#include <tuple>
#include <string>

#include <algorithm>
#include <numeric>

#include <iostream>
#include <fstream>
#include <iomanip>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include "maths.h"
#include "units.h"
#include "phys.h"
#include "algos.h"
#include "expr.h"


namespace tl2_mag {

// ----------------------------------------------------------------------------
// helper functions
// ----------------------------------------------------------------------------

/**
 * rotate spin vector for incommensurate structures, i.e. helices
 */
template<class t_mat, class t_vec, class t_real = typename t_mat::value_type>
requires tl2::is_mat<t_mat> && tl2::is_vec<t_vec>
void rotate_spin_incommensurate(t_vec& spin_vec,
	const t_vec& sc_vec, const t_vec& ordering, const t_vec& rotaxis,
	t_real eps = std::numeric_limits<t_real>::epsilon())
{
	t_real sc_angle = t_real(2) * tl2::pi<t_real> *
		tl2::inner<t_vec>(ordering, sc_vec);

	if(!tl2::equals_0<t_real>(sc_angle, t_real(eps)))
	{
		t_mat sc_rot = tl2::rotation<t_mat, t_vec>(rotaxis, sc_angle);
		spin_vec = sc_rot * spin_vec;
	}
}


/**
 * polarisation matrix
 * @see https://doi.org/10.1088/1361-6463/aa7573
 */
template<class t_mat, class t_cplx = typename t_mat::value_type>
requires tl2::is_mat<t_mat>
t_mat get_polarisation(int channel = 0, bool in_chiral_base = true)
{
	if(in_chiral_base)
	{
		t_mat pol = tl2::zero<t_mat>(3);

		// just pick the selected component on the diagonal
		if(channel >=0 && channel < 3)
			pol(channel, channel) = t_cplx(1);

		return pol;
	}
	else
	{
		constexpr const t_cplx halfi = t_cplx(0, 0.5);
		constexpr const t_cplx half = t_cplx(0.5, 0);

		// TODO: check coordinate system
		switch(channel)
		{
			case 0: return tl2::create<t_mat>({
				{   half, +halfi, 0 },
				{ -halfi,   half, 0 },
				{      0,      0, 0 } });
			case 1: return tl2::create<t_mat>({
				{   half, -halfi, 0 },
				{ +halfi,   half, 0 },
				{      0,      0, 0 } });
			case 2: return tl2::create<t_mat>({
				{      0,      0, 0 },
				{      0,      0, 0 },
				{      0,      0, 1 } });
		}
	}

	return tl2::zero<t_mat>(3);
}
// ----------------------------------------------------------------------------



/**
 * calculates magnon dynamics,
 * implementing the formalism given in (Toth 2015)
 */
template<
	class t_mat, class t_vec,
	class t_mat_real, class t_vec_real,
	class t_cplx = typename t_mat::value_type,
	class t_real = typename t_mat_real::value_type,
	class t_size = std::size_t>
requires tl2::is_mat<t_mat> && tl2::is_vec<t_vec> &&
	tl2::is_mat<t_mat_real> && tl2::is_vec<t_vec_real>
class MagDyn
{
public:
	// ----------------------------------------------------------------------------
	// input- and output structs
	// ----------------------------------------------------------------------------
	/**
	 * magnetic atom sites
	 */
	struct AtomSite
	{
		std::string name{};      // identifier
		t_size index{};          // index

		t_vec_real pos{};        // atom position

		std::string spin_dir[3]; // expression for spin direction
		t_real spin_mag{};       // spin magnitude
		t_mat g{};               // g factor
	};


	/**
	 * temporary per-site calculation results
	 */
	struct AtomSiteCalc
	{
		t_vec spin_dir{};        // spin direction

		t_vec u{}, u_conj{};
		t_vec v{};
	};


	/**
	 * couplings between magnetic atoms
	 */
	struct ExchangeTerm
	{
		std::string name{};      // identifier
		t_size index{};          // index

		t_size atom1{};          // atom 1 index
		t_size atom2{};          // atom 2 index
		t_vec_real dist{};       // distance between unit cells

		std::string J{};         // parsable expression for Heisenberg interaction
		std::string dmi[3];      // parsable expression for Dzyaloshinskij-Moriya interaction
	};


	/**
	 * temporary per-term calculation results
	 */
	struct ExchangeTermCalc
	{
		t_cplx J{};              // Heisenberg interaction
		t_vec dmi{};             // Dzyaloshinskij-Moriya interaction
	};


	/**
	 * terms related to an external magnetic field
	 */
	struct ExternalField
	{
		bool align_spins{};      // align spins along external field
		t_vec_real dir{};        // field direction
		t_real mag{};            // field magnitude
	};


	/**
	 * eigenenergies and spin-spin correlation matrix
	 */
	struct EnergyAndWeight
	{
		t_real E{};

		// full dynamical structure factor
		t_mat S{};
		t_real weight_full{};
		t_real weight_channel_full[3] = {0., 0., 0.};

		// projected dynamical structure factor for neutron scattering
		t_mat S_perp{};
		t_real weight{};
		t_real weight_channel[3] = {0., 0., 0.};
	};


	/**
	 * variables for the expression parser
	 */
	struct Variable
	{
		std::string name{};
		t_cplx value{};
	};
	// ----------------------------------------------------------------------------


public:
	MagDyn() = default;
	~MagDyn() = default;

	MagDyn(const MagDyn&) = default;
	MagDyn& operator=(const MagDyn&) = default;


	// --------------------------------------------------------------------
	// cleanup functions
	// --------------------------------------------------------------------
	/**
	 * clear all
	 */
	void Clear()
	{
		ClearVariables();
		ClearAtomSites();
		ClearExchangeTerms();
		ClearExternalField();

		// clear temperature, -1: don't use
		m_temperature = -1.;

		// clear ordering wave vector
		m_ordering = tl2::zero<t_vec_real>(3);

		// reset rotation axis
		m_rotaxis = tl2::create<t_vec_real>({1., 0., 0.});
	}


	/**
	 * clear all parser variables
	 */
	void ClearVariables()
	{
		m_variables.clear();
	}


	/**
	 * clear all atom sites
	 */
	void ClearAtomSites()
	{
		m_sites.clear();
		m_sites_calc.clear();
	}


	/**
	 * clear all couplings
	 */
	void ClearExchangeTerms()
	{
		m_exchange_terms.clear();
		m_exchange_terms_calc.clear();
	}


	/**
	 * clear the external field settings
	 */
	void ClearExternalField()
	{
		m_field.dir.clear();
		m_field.mag = 0.;
		m_field.align_spins = false;
	}
	// --------------------------------------------------------------------


	// --------------------------------------------------------------------
	// getter
	// --------------------------------------------------------------------
	const std::vector<Variable>& GetVariables() const { return m_variables; }
	const std::vector<AtomSite>& GetAtomSites() const { return m_sites; }
	const std::vector<AtomSiteCalc>& GetAtomSitesCalc() const { return m_sites_calc; }
	const std::vector<ExchangeTerm>& GetExchangeTerms() const { return m_exchange_terms; }
	const std::vector<ExchangeTermCalc>& GetExchangeTermsCalc() const { return m_exchange_terms_calc; }

	const ExternalField& GetExternalField() const { return m_field; }
	const t_vec_real& GetRotationAxis() const { return m_rotaxis; }
	const t_vec_real& GetOrderingWavevector() const { return m_ordering; }

	t_real GetTemperature() const { return m_temperature; }
	t_real GetBoseCutoffEnergy() const { return m_bose_cutoff; }

	bool IsIncommensurate() const
	{
		return m_is_incommensurate || m_force_incommensurate;
	}
	// --------------------------------------------------------------------


	// --------------------------------------------------------------------
	// setter
	// --------------------------------------------------------------------
	void SetEpsilon(t_real eps) { m_eps = eps; }
	void SetPrecision(int prec) { m_prec = prec; }

		void SetTemperature(t_real T) { m_temperature = T; }
	void SetBoseCutoffEnergy(t_real E) { m_bose_cutoff = E; }
	void SetUniteDegenerateEnergies(bool b) { m_unite_degenerate_energies = b; }
	void SetForceIncommensurate(bool b) { m_force_incommensurate = b; }


	void SetExternalField(const ExternalField& field)
	{
		m_field = field;
		t_real len = tl2::norm<t_vec_real>(m_field.dir);
		if(!tl2::equals_0<t_real>(len, m_eps))
			m_field.dir /= len;
	}


	void SetOrderingWavevector(const t_vec_real& ordering)
	{
		m_ordering = ordering;
		m_is_incommensurate = !tl2::equals_0<t_vec_real>(m_ordering, m_eps);
	}


	void SetRotationAxis(const t_vec_real& axis)
	{
		m_rotaxis = axis;
		t_real len = tl2::norm<t_vec_real>(m_rotaxis);
		if(!tl2::equals_0<t_real>(len, m_eps))
			m_rotaxis /= len;
	}


	void AddVariable(Variable&& var)
	{
		m_variables.emplace_back(std::forward<Variable&&>(var));
	}


	void SetVariable(Variable&& var)
	{
		// is a variable with the same name already registered?
		auto iter = std::find_if(m_variables.begin(), m_variables.end(),
			[&var](const auto& thevar)
		{
			return thevar.name == var.name;
		});

		if(iter == m_variables.end())
		{
			// add a new variable
			AddVariable(std::forward<Variable&&>(var));
		}
		else
		{
			// replace the value of an existing variable
			iter->value = var.value;
		}
	}


	void AddAtomSite(AtomSite&& site, bool set_index = true)
	{
		if(set_index)
			site.index = GetAtomSites().size();
		m_sites.emplace_back(std::forward<AtomSite&&>(site));
	}


	void AddExchangeTerm(ExchangeTerm&& term, bool set_index = true)
	{
		if(set_index)
			term.index = GetExchangeTerms().size();
		m_exchange_terms.emplace_back(std::forward<ExchangeTerm&&>(term));
	}
	// --------------------------------------------------------------------


	// --------------------------------------------------------------------
	/**
	 * get an expression parser object with registered variables
	 */
	tl2::ExprParser<t_cplx> GetExprParser() const
	{
		tl2::ExprParser<t_cplx> parser;

		// register all variables
		parser.SetAutoregisterVariables(false);
		for(const Variable& var : m_variables)
			parser.register_var(var.name, var.value);

		return parser;
	}
	// --------------------------------------------------------------------


	/**
	 * get the needed supercell ranges from the exchange terms
	 */
	std::tuple<t_vec_real, t_vec_real> GetSupercellMinMax() const
	{
		t_vec_real min = tl2::zero<t_vec_real>(3);
		t_vec_real max = tl2::zero<t_vec_real>(3);

		for(const ExchangeTerm& term : m_exchange_terms)
		{
			for(t_size i=0; i<3; ++i)
			{
				min[i] = std::min(min[i], term.dist[i]);
				max[i] = std::max(max[i], term.dist[i]);
			}
		}

		return std::make_tuple(min, max);
	}


	/**
	 * calculate the spin rotation trafo for the atom sites
	 */
	void CalcAtomSites()
	{
		const t_size num_sites = m_sites.size();
		if(num_sites == 0)
			return;

		m_sites_calc.clear();
		m_sites_calc.reserve(num_sites);

		bool use_field =
			(!tl2::equals_0<t_real>(m_field.mag, m_eps) || m_field.align_spins)
			&& m_field.dir.size() == 3;

		if(use_field)
		{
			// rotate field to [001] direction
			m_rot_field = tl2::convert<t_mat>(
				tl2::trans<t_mat_real>(
				tl2::rotation<t_mat_real, t_vec_real>(
				-m_field.dir, m_zdir, &m_rotaxis, m_eps)));

			/*std::cout << "Field rotation from:\n";
			tl2::niceprint(std::cout, -m_field.dir, 1e-4, 4);
			std::cout << "\nto:\n";
			tl2::niceprint(std::cout, m_zdir, 1e-4, 4);
			std::cout << "\nmatrix:\n";
			tl2::niceprint(std::cout, m_rot_field, 1e-4, 4);
			std::cout << std::endl;*/
		}

		tl2::ExprParser parser = GetExprParser();

		for(t_size site_idx=0; site_idx<m_sites.size(); ++site_idx)
		{
			try
			{
				const AtomSite& site = m_sites[site_idx];
				AtomSiteCalc site_calc{};

				site_calc.spin_dir = tl2::zero<t_vec>(3);

				for(t_size dir_idx=0; dir_idx<3; ++dir_idx)
				{
					// empty string?
					if(!site.spin_dir[dir_idx].size())
						continue;

					if(bool dir_ok = parser.parse(site.spin_dir[dir_idx]); dir_ok)
					{
						site_calc.spin_dir[dir_idx] = parser.eval();
					}
					else
					{
						std::cerr << "Error parsing in spin direction \""
							<< site.spin_dir[dir_idx]
							<< "\"." << std::endl;
					}
				}

				//tl2::niceprint(std::cout, site_calc.spin_dir, 1e-4, 4);
				//std::cout << std::endl;

				// spin rotation of equation (9) from (Toth 2015)
				if(m_field.align_spins)
				{
					std::tie(site_calc.u, site_calc.v) =
						R_to_uv(m_rot_field);
				}
				else
				{
					std::tie(site_calc.u, site_calc.v) =
						spin_to_uv(site_calc.spin_dir);
				}

				site_calc.u_conj = tl2::conj(site_calc.u);
				m_sites_calc.emplace_back(std::move(site_calc));
			}
			catch(const std::exception& ex)
			{
				std::cerr << ex.what() << std::endl;
			}
		}
	}


	/**
	 * parse the exchange term expressions
	 */
	void CalcExchangeTerms()
	{
		const t_size num_terms = m_exchange_terms.size();
		if(num_terms == 0)
			return;

		tl2::ExprParser parser = GetExprParser();

		m_exchange_terms_calc.clear();
		m_exchange_terms_calc.reserve(num_terms);

		for(t_size term_idx=0; term_idx<num_terms; ++term_idx)
		{
			const ExchangeTerm& term = m_exchange_terms[term_idx];
			ExchangeTermCalc calc;

			try
			{
				if(bool J_ok = parser.parse(term.J); J_ok)
				{
					t_cplx J = parser.eval();
					calc.J = J;
				}
				else
				{
					std::cerr << "Error parsing J term \""
						<< term.J << "\"."
						<< std::endl;
				}

				calc.dmi = tl2::zero<t_vec>(3);

				for(t_size dmi_idx=0; dmi_idx<3; ++dmi_idx)
				{
					// empty string?
					if(!term.dmi[dmi_idx].size())
						continue;

					if(bool dmi_ok = parser.parse(term.dmi[dmi_idx]); dmi_ok)
					{
						calc.dmi[dmi_idx] = parser.eval();
					}
					else
					{
						std::cerr << "Error parsing DMI term \""
							<< term.dmi[dmi_idx]
							<< "\"." << std::endl;
					}
				}
			}
			catch(const std::exception& ex)
			{
				std::cerr << ex.what() << std::endl;
			}

			m_exchange_terms_calc.emplace_back(std::move(calc));
		}
	}


	/**
	 * update the site and term indices
	 */
	void CalcIndices()
	{
		for(t_size site_idx=0; site_idx<m_sites.size(); ++site_idx)
		{
			AtomSite& site = m_sites[site_idx];
			site.index = site_idx;
		}

		for(t_size term_idx=0; term_idx<m_exchange_terms.size(); ++term_idx)
		{
			ExchangeTerm& term = m_exchange_terms[term_idx];
			term.index = term_idx;
		}
	}


	/**
	 * get the hamiltonian at the given momentum
	 * (CalcAtomSites() needs to be called once before this function)
	 * @note implements the formalism given by (Toth 2015)
	 * @note a first version for a simplified ferromagnetic dispersion was based on (Heinsdorf 2021)
	 */
	t_mat GetHamiltonian(const t_vec_real& Qvec) const
	{
		const t_size num_sites = m_sites.size();
		const t_size num_terms = m_exchange_terms.size();

		if(num_sites == 0 || num_terms == 0)
			return {};
		if(num_sites != m_sites_calc.size() ||
			num_terms != m_exchange_terms_calc.size())
			return {};

		// bohr magneton in [meV/T]
		constexpr const t_real muB = tl2::muB<t_real>
			/ tl2::meV<t_real> * tl2::tesla<t_real>;

		// build the interaction matrices J(Q) and J(-Q) of
		// equations (12) and (14) from (Toth 2015)
		t_mat J_Q = tl2::zero<t_mat>(num_sites*3, num_sites*3);
		t_mat J_mQ = tl2::zero<t_mat>(num_sites*3, num_sites*3);
		t_mat J_Q0 = tl2::zero<t_mat>(num_sites*3, num_sites*3);

		// iterate couplings
		for(t_size term_idx=0; term_idx<num_terms; ++term_idx)
		{
			const ExchangeTerm& term = m_exchange_terms[term_idx];
			const ExchangeTermCalc& term_calc = m_exchange_terms_calc[term_idx];

			if(term.atom1 >= num_sites || term.atom2 >= num_sites)
			{
				std::cerr << "Error: Site index out of bounds for coupling term "
					<< term_idx << "." << std::endl;
				continue;
			}

			// exchange interaction matrix with dmi as anti-symmetric part, see (Toth 2015) p. 2
			t_mat J = tl2::diag<t_mat>(
				tl2::create<t_vec>({term_calc.J, term_calc.J, term_calc.J}));

			if(term_calc.dmi.size() == 3)
			{
				// cross product matrix
				J += tl2::skewsymmetric<t_mat, t_vec>(-term_calc.dmi);
			}

			// incommensurate case: rotation wrt magnetic unit cell
			// equations (21), (6), (2) as well as section 10 from (Toth 2015)
			if(IsIncommensurate())
			{
				t_real rot_UC_angle = s_twopi * tl2::inner<t_vec_real>(m_ordering, term.dist);
				if(!tl2::equals_0<t_real>(rot_UC_angle, m_eps))
				{
					t_mat rot_UC = tl2::convert<t_mat>(
						 tl2::rotation<t_mat_real, t_vec_real>(
							m_rotaxis, rot_UC_angle));
					J = J * rot_UC;
				}
			}

			// equation (14) from (Toth 2015)
			t_real phase_sign = -1.;
			t_cplx phase_Q = std::exp(phase_sign * s_imag * s_twopi *
				tl2::inner<t_vec_real>(term.dist, Qvec));
			t_cplx phase_mQ = std::exp(phase_sign * s_imag * s_twopi *
				tl2::inner<t_vec_real>(-term.dist, Qvec));

			t_mat J_T = tl2::trans(J);
			t_real factor = /*0.5*/ 1.;

			// include these two terms to fulfill equation (11) from (Toth 2015)
			tl2::add_submat<t_mat>(J_Q, factor * J * phase_Q, term.atom1*3, term.atom2*3);
			tl2::add_submat<t_mat>(J_Q, factor * J_T * phase_mQ, term.atom2*3, term.atom1*3);

			tl2::add_submat<t_mat>(J_mQ, factor * J * phase_mQ, term.atom1*3, term.atom2*3);
			tl2::add_submat<t_mat>(J_mQ, factor * J_T * phase_Q, term.atom2*3, term.atom1*3);

			tl2::add_submat<t_mat>(J_Q0, factor * J, term.atom1*3, term.atom2*3);
			tl2::add_submat<t_mat>(J_Q0, factor * J_T, term.atom2*3, term.atom1*3);
		}


		// create the hamiltonian of equation (25) and (26) from (Toth 2015)
		t_mat A = tl2::create<t_mat>(num_sites, num_sites);
		t_mat A_mQ = tl2::create<t_mat>(num_sites, num_sites);
		t_mat B = tl2::create<t_mat>(num_sites, num_sites);
		t_mat C = tl2::zero<t_mat>(num_sites, num_sites);

		bool use_field = !tl2::equals_0<t_real>(m_field.mag, m_eps)
			&& m_field.dir.size() == 3;

		// iterate sites
		for(t_size i=0; i<num_sites; ++i)
		for(t_size j=0; j<num_sites; ++j)
		{
			t_mat J_sub_Q = tl2::submat<t_mat>(J_Q, i*3, j*3, 3, 3);
			t_mat J_sub_mQ = tl2::submat<t_mat>(J_mQ, i*3, j*3, 3, 3);

			t_real S_i = m_sites[i].spin_mag;
			t_real S_j = m_sites[j].spin_mag;

			// get the precalculated u and v vectors for the commensurate case
			const t_vec& u_i = m_sites_calc[i].u;
			const t_vec& u_j = m_sites_calc[j].u;
			const t_vec& u_conj_j = m_sites_calc[j].u_conj;
			const t_vec& v_i = m_sites_calc[i].v;

			// equation (26) from (Toth 2015)
			t_real SiSj = 0.5 * std::sqrt(S_i*S_j);
			A(i, j) = SiSj * tl2::inner_noconj<t_vec>(u_i, J_sub_Q * u_conj_j);
			A_mQ(i, j) = SiSj * tl2::inner_noconj<t_vec>(u_i, J_sub_mQ * u_conj_j);
			B(i, j) = SiSj * tl2::inner_noconj<t_vec>(u_i, J_sub_Q * u_j);

			if(i == j)
			{
				for(t_size k=0; k<num_sites; ++k)
				{
					t_real S_k = m_sites[k].spin_mag;

					// get the precalculated v_k vectors for the commensurate case
					const t_vec& v_k = m_sites_calc[k].v;

					// equation (26) from (Toth 2015)
					t_mat J_sub_Q0 = tl2::submat<t_mat>(J_Q0, i*3, k*3, 3, 3);
					C(i, j) += S_k * tl2::inner_noconj<t_vec>(v_i, J_sub_Q0 * v_k);
				}

				// include external field, equation (28) from (Toth 2015)
				if(use_field)
				{
					t_vec B = tl2::convert<t_vec>(-m_field.dir) * m_field.mag;

					t_vec gv = m_sites[i].g * v_i;
					t_cplx Bgv = tl2::inner_noconj<t_vec>(B, gv);

					A(i, j) -= 0.5 * muB * Bgv;
					A_mQ(i, j) -= 0.5 * muB * Bgv;
				}
			}
		}

		// equation (25) from (Toth 2015)
		t_mat H = tl2::zero<t_mat>(num_sites*2, num_sites*2);
		tl2::set_submat(H, A - C, 0, 0);
		tl2::set_submat(H, B, 0, num_sites);
		tl2::set_submat(H, tl2::herm(B), num_sites, 0);
		tl2::set_submat(H, tl2::conj(A_mQ) - C, num_sites, num_sites);

		return H;
	}


	/**
	 * get the energies and the dynamical structure factor from a hamiltonian
	 * @note implements the formalism given by (Toth 2015)
	 */
	std::vector<EnergyAndWeight> GetEnergiesFromHamiltonian(t_mat _H, const t_vec_real& Qvec,
		bool only_energies = false) const
	{
		const t_size num_sites = m_sites.size();
		if(num_sites == 0 || _H.size1() == 0)
			return {};

		// equation (30) from (Toth 2015)
		t_mat g_sign = tl2::zero<t_mat>(num_sites*2, num_sites*2);
		for(t_size i=0; i<num_sites; ++i)
			g_sign(i, i) = 1.;
		for(t_size i=num_sites; i<2*num_sites; ++i)
			g_sign(i, i) = -1.;

		// equation (31) from (Toth 2015)
		t_mat C_mat;
		t_size chol_try = 0;
		for(; chol_try<m_tries_chol; ++chol_try)
		{
			auto [chol_ok, _C] = tl2_la::chol<t_mat>(_H);

			if(chol_ok)
			{
				C_mat = std::move(_C);
				break;
			}
			else
			{
				if(chol_try == m_tries_chol-1)
				{
					using namespace tl2_ops;
					std::cerr << "Warning: Cholesky decomposition failed for Q = "
						<< Qvec << "." << std::endl;
					C_mat = std::move(_C);
					break;
				}

				// try forcing the hamilton to be positive definite
				for(t_size i=0; i<2*num_sites; ++i)
					_H(i, i) += m_delta_chol;
			}
		}

		if(chol_try > 1)
		{
			using namespace tl2_ops;
			std::cerr << "Warning: Needed " << chol_try
				<< " corrections for cholesky decomposition for Q = "
				<< Qvec << "." << std::endl;
		}

		t_mat C_herm = tl2::herm<t_mat>(C_mat);

		// see p. 5 in (Toth 2015)
		t_mat H_mat = C_mat * g_sign * C_herm;
		//tl2::niceprint(std::cout, H, 1e-4, 4);
		//std::cout << std::endl;

		bool is_herm = tl2::is_symm_or_herm<t_mat, t_real>(H_mat, m_eps);
		if(!is_herm)
		{
			using namespace tl2_ops;
			std::cerr << "Warning: Hamiltonian is not hermitian for Q = "
				<< Qvec << "." << std::endl;
		}

		// eigenvalues of the hamiltonian correspond to the energies
		// eigenvectors correspond to the spectral weights
		auto [evecs_ok, evals, evecs] =
			tl2_la::eigenvec<t_mat, t_vec, t_cplx, t_real>(
				H_mat, only_energies, is_herm, true);
		if(!evecs_ok)
		{
			using namespace tl2_ops;
			std::cerr << "Warning: Eigensystem calculation failed for Q = "
				<< Qvec << "." << std::endl;
		}


		std::vector<EnergyAndWeight> energies_and_correlations{};
		energies_and_correlations.reserve(evals.size());

		// register energies
		for(const auto& eval : evals)
		{
			EnergyAndWeight EandS { .E = eval.real(), };
			energies_and_correlations.emplace_back(std::move(EandS));
		}


		// weight factors
		if(!only_energies)
		{
			// get the sorting of the energies
			std::vector<t_size> sorting = tl2::get_perm(
				energies_and_correlations.size(),
				[&energies_and_correlations](t_size idx1, t_size idx2) -> bool
				{
					return energies_and_correlations[idx1].E >=
						energies_and_correlations[idx2].E;
				});

			//energies_and_correlations = tl2::reorder(energies_and_correlations, sorting);
			evecs = tl2::reorder(evecs, sorting);
			evals = tl2::reorder(evals, sorting);

			/*for(std::size_t idx=0; idx<evecs.size(); ++idx)
			{
				std::cout << "eval = " << evals[idx] << std::endl;
				tl2::niceprint(std::cout, evecs[idx], 1e-4, 4);
				std::cout << std::endl;
			}*/

			t_mat evec_mat = tl2::create<t_mat>(evecs);
			t_mat evec_mat_herm = tl2::herm(evec_mat);
			//tl2::niceprint(std::cout, evec_mag, 1e-4, 4);
			//std::cout << std::endl;

			// equation (32) from (Toth 2015)
			t_mat L_mat = evec_mat_herm * H_mat * evec_mat; // energies
			t_mat E_sqrt = g_sign * L_mat;                  // abs. energies
			for(t_size i=0; i<E_sqrt.size1(); ++i)
				E_sqrt(i, i) = std::sqrt(E_sqrt/*L_mat*/(i, i)); // sqrt. of abs. energies
			//tl2::niceprint(std::cout, L, 1e-4, 4);
			//std::cout << std::endl;

			// re-create energies, to be consistent with the weights
			energies_and_correlations.clear();
			for(t_size i=0; i<L_mat.size1(); ++i)
			{
				EnergyAndWeight EandS
				{
					.E = L_mat(i, i).real(),
					.S = tl2::zero<t_mat>(3, 3),
					.S_perp = tl2::zero<t_mat>(3, 3),
				};

				energies_and_correlations.emplace_back(std::move(EandS));
			}

			auto [C_inv, inv_ok] = tl2::inv(C_mat);
			if(!inv_ok)
			{
				using namespace tl2_ops;
				std::cerr << "Warning: Inversion failed for Q = "
					<< Qvec << "." << std::endl;
			}

			// equation (34) from (Toth 2015)
			t_mat trafo = C_inv * evec_mat * E_sqrt;
			t_mat trafo_herm = tl2::herm(trafo);

			//t_mat D = trafo_herm * _H * trafo;
			//tl2::niceprint(std::cout, D, 1e-4, 4);
			//tl2::niceprint(std::cout, E, 1e-4, 4);
			//tl2::niceprint(std::cout, L, 1e-4, 4);
			//std::cout << std::endl;

			/*std::cout << "Y = np.zeros(3*3*4*4, dtype=complex).reshape((4,4,3,3))" << std::endl;
			std::cout << "V = np.zeros(3*3*4*4, dtype=complex).reshape((4,4,3,3))" << std::endl;
			std::cout << "Z = np.zeros(3*3*4*4, dtype=complex).reshape((4,4,3,3))" << std::endl;
			std::cout << "W = np.zeros(3*3*4*4, dtype=complex).reshape((4,4,3,3))" << std::endl;*/

			// building the spin correlation functions of equation (47) from (Toth 2015)
			for(int x_idx=0; x_idx<3; ++x_idx)
			for(int y_idx=0; y_idx<3; ++y_idx)
			{
				// equations (44) from (Toth 2015)
				auto create_matrices = [num_sites](
					t_mat& V, t_mat& W, t_mat& Y, t_mat& Z)
				{
					V = tl2::create<t_mat>(num_sites, num_sites);
					W = tl2::create<t_mat>(num_sites, num_sites);
					Y = tl2::create<t_mat>(num_sites, num_sites);
					Z = tl2::create<t_mat>(num_sites, num_sites);
				};

				t_mat V, W, Y, Z;
				create_matrices(V, W, Y, Z);

				for(t_size i=0; i<num_sites; ++i)
				for(t_size j=0; j<num_sites; ++j)
				{
					auto calc_mat_elems = [this, i, j, x_idx, y_idx]
						(const t_vec_real& Qvec,
						 t_mat& Y, t_mat& V, t_mat& Z, t_mat& W)
					{
						const t_vec_real& pos_i = m_sites[i].pos;
						const t_vec_real& pos_j = m_sites[j].pos;

						t_real S_i = m_sites[i].spin_mag;
						t_real S_j = m_sites[j].spin_mag;

						// get the precalculated u vectors
						const t_vec& u_i = m_sites_calc[i].u;
						const t_vec& u_j = m_sites_calc[j].u;
						const t_vec& u_conj_i = m_sites_calc[i].u_conj;
						const t_vec& u_conj_j = m_sites_calc[j].u_conj;

						// TODO: check these
						t_real SiSj = 4. * std::sqrt(S_i*S_j);
						t_real phase_sign = 1.;

						t_cplx phase = std::exp(phase_sign * s_imag * s_twopi *
							tl2::inner<t_vec_real>(pos_j - pos_i, Qvec));
						phase *= SiSj;

						// matrix elements of equ. (44) from (Toth 2015)
						Y(i, j) = phase * u_i[x_idx] * u_conj_j[y_idx];
						V(i, j) = phase * u_conj_i[x_idx] * u_conj_j[y_idx];
						Z(i, j) = phase * u_i[x_idx] * u_j[y_idx];
						W(i, j) = phase * u_conj_i[x_idx] * u_j[y_idx];

						/*std::cout << "Y[" << i << ", " << j << ", "
							<< x_idx << ", " << y_idx << "] = "
							<< Y(i, j).real() << " + " << Y(i, j).imag() << "j\n"
							<< "V[" << i << ", " << j << ", "
							<< x_idx << ", " << y_idx << "] = "
							<< V(i, j).real() << " + " << V(i, j).imag() << "j\n"
							<< "Z[" << i << ", " << j << ", "
							<< x_idx << ", " << y_idx << "] = "
							<< Z(i, j).real() << " + " << Z(i, j).imag() << "j\n"
							<< "W[" << i << ", " << j << ", "
							<< x_idx << ", " << y_idx << "] = "
							<< W(i, j).real() << " + " << W(i, j).imag() << "j"
							<< std::endl;*/
					};

					calc_mat_elems(Qvec, Y, V, Z, W);
				} // end of iteration over sites


				auto calc_S = [this, num_sites, x_idx, y_idx, &trafo, &trafo_herm, &energies_and_correlations]
					(t_mat EnergyAndWeight::*S, const t_mat& Y, const t_mat& V, const t_mat& Z, const t_mat& W)
				{
					// equation (47) from (Toth 2015)
					t_mat M = tl2::create<t_mat>(num_sites*2, num_sites*2);
					tl2::set_submat(M, Y, 0, 0);
					tl2::set_submat(M, V, num_sites, 0);
					tl2::set_submat(M, Z, 0, num_sites);
					tl2::set_submat(M, W, num_sites, num_sites);

					t_mat M_trafo = trafo_herm * M * trafo;

					/*std::cout << "M_trafo for x=" << x_idx << ", y=" << y_idx << ":\n";
					tl2::niceprint(std::cout, M_trafo, 1e-4, 4);
					std::cout << std::endl;*/

					for(t_size i=0; i<energies_and_correlations.size(); ++i)
						(energies_and_correlations[i].*S)(x_idx, y_idx) += M_trafo(i, i) / t_real(2*num_sites);
				};

				calc_S(&EnergyAndWeight::S, Y, V, Z, W);
			} // end of coordinate iteration
		} // end of weight calculation

		return energies_and_correlations;
	}


	/**
	 * applies projectors and weight factors to get neutron intensities
	 * @note implements the formalism given by (Toth 2015)
	 */
	void GetIntensities(const t_vec_real& Qvec, std::vector<EnergyAndWeight>& energies_and_correlations) const
	{
		for(EnergyAndWeight& E_and_S : energies_and_correlations)
		{
			// apply bose factor
			if(m_temperature >= 0.)
				E_and_S.S *= tl2::bose_cutoff(E_and_S.E, m_temperature, m_bose_cutoff);

			// apply orthogonal projector for magnetic neutron scattering,
			// see (Shirane 2002), p. 37, eq. (2.64)
			//t_vec bragg_rot = use_field ? m_rot_field * m_bragg : m_bragg;
			//proj_neutron = tl2::ortho_projector<t_mat, t_vec>(bragg_rot, false);
			t_mat proj_neutron = tl2::ortho_projector<t_mat, t_vec>(Qvec, false);
			E_and_S.S_perp = proj_neutron * E_and_S.S * proj_neutron;

			// weights
			E_and_S.weight = std::abs(tl2::trace<t_mat>(E_and_S.S_perp).real());
			E_and_S.weight_full = std::abs(tl2::trace<t_mat>(E_and_S.S).real());

			// polarisation channels
			for(int i=0; i<3; ++i)
			{
				t_mat pol = get_polarisation<t_mat>(i, false);
				t_mat Sperp = pol * E_and_S.S_perp;
				t_mat S = pol * E_and_S.S;

				E_and_S.weight_channel[i] = std::abs(tl2::trace<t_mat>(Sperp).real());
				E_and_S.weight_channel_full[i] = std::abs(tl2::trace<t_mat>(S).real());
			}
		}
	}


	/**
	 * unite degenerate energies and their corresponding eigenstates
	 */
	std::vector<EnergyAndWeight> UniteEnergies(const std::vector<EnergyAndWeight>& energies_and_correlations) const
	{
		std::vector<EnergyAndWeight> new_energies_and_correlations{};
		new_energies_and_correlations.reserve(energies_and_correlations.size());

		for(const auto& curState : energies_and_correlations)
		{
			t_real curE = curState.E;

			auto iter = std::find_if(
				new_energies_and_correlations.begin(),
					new_energies_and_correlations.end(),
					[curE, this](const auto& E_and_S) -> bool
					{
						t_real E = E_and_S.E;
						return tl2::equals<t_real>(E, curE, m_eps);
					});

			if(iter == new_energies_and_correlations.end())
			{
				// energy not yet seen
				new_energies_and_correlations.push_back(curState);
			}
			else
			{
				// energy already seen: add correlation matrices and weights
				iter->S += curState.S;
				iter->S_perp += curState.S_perp;

				iter->weight += curState.weight;
				iter->weight_full += curState.weight_full;

				for(int i=0; i<3; ++i)
				{
					iter->weight_channel[i] += curState.weight_channel[i];
					iter->weight_channel_full[i] += curState.weight_channel_full[i];
				}
			}
		}

		return new_energies_and_correlations;
	}


	/**
	 * get the energies and the spin-correlation at the given momentum
	 * (also calculates incommensurate contributions and applies weight factors)
	 * @note implements the formalism given by (Toth 2015)
	 */
	std::vector<EnergyAndWeight> GetEnergies(const t_vec_real& Qvec,
		bool only_energies = false) const
	{
		t_mat H = GetHamiltonian(Qvec);
		std::vector<EnergyAndWeight> EandWs = GetEnergiesFromHamiltonian(H, Qvec, only_energies);

		if(IsIncommensurate())
		{
			// equations 39 and 40 from (Toth 2015)
			t_mat proj_norm = tl2::convert<t_mat>(
				tl2::projector<t_mat_real, t_vec_real>(
					m_rotaxis, true));

			t_real phase_sign = -1.;
			t_mat rot_incomm = tl2::unit<t_mat>(3);
			rot_incomm -= s_imag * phase_sign * tl2::skewsymmetric<t_mat, t_vec>(m_rotaxis);
			rot_incomm -= proj_norm;
			rot_incomm *= 0.5;

			t_mat rot_incomm_conj = tl2::conj(rot_incomm);

			t_mat H_p = GetHamiltonian(Qvec + m_ordering);
			t_mat H_m = GetHamiltonian(Qvec - m_ordering);

			std::vector<EnergyAndWeight> EandWs_p = GetEnergiesFromHamiltonian(
				H_p, Qvec + m_ordering, only_energies);
			std::vector<EnergyAndWeight> EandWs_m = GetEnergiesFromHamiltonian(
				H_m, Qvec - m_ordering, only_energies);

			if(!only_energies)
			{
				// formula 40 from (Toth 2015)
				for(EnergyAndWeight& EandW : EandWs)
					EandW.S = EandW.S * proj_norm;
				for(EnergyAndWeight& EandW : EandWs_p)
					EandW.S = EandW.S * rot_incomm;
				for(EnergyAndWeight& EandW : EandWs_m)
					EandW.S = EandW.S * rot_incomm_conj;
			}

			// unite energies and weights
			for(EnergyAndWeight& EandW : EandWs_p)
				EandWs.emplace_back(std::move(EandW));
			for(EnergyAndWeight& EandW : EandWs_m)
				EandWs.emplace_back(std::move(EandW));
		}

		if(!only_energies)
			GetIntensities(Qvec, EandWs);

		if(m_unite_degenerate_energies)
			EandWs = UniteEnergies(EandWs);

		return EandWs;
	}


	std::vector<EnergyAndWeight> GetEnergies(t_real h, t_real k, t_real l,
		bool only_energies = false) const
	{
		// momentum transfer
		const t_vec_real Qvec = tl2::create<t_vec_real>({ h, k, l });
		return GetEnergies(Qvec, only_energies);
	}


	/**
	 * get the energy of the goldstone mode
	 * @note a first version for a simplified ferromagnetic dispersion was based on (Heinsdorf 2021)
	 */
	t_real GetGoldstoneEnergy() const
	{
		auto energies_and_correlations = GetEnergies(0., 0., 0., true);
		auto min_iter = std::min_element(
			energies_and_correlations.begin(), energies_and_correlations.end(),
			[](const auto& E_and_S_1, const auto& E_and_S_2) -> bool
			{
				return E_and_S_1.E < E_and_S_2.E;
			});

		if(min_iter != energies_and_correlations.end())
			return min_iter->E;

		return 0.;
	}


	/**
	 * generates the dispersion plot along the given q path
	 */
	void SaveDispersion(const std::string& filename,
		t_real h_start, t_real k_start, t_real l_start,
		t_real h_end, t_real k_end, t_real l_end,
		t_size num_qs = 128) const
	{
		std::ofstream ofstr{filename};
		ofstr.precision(m_prec);

		ofstr
			<< std::setw(m_prec*2) << std::left << "# h"
			<< std::setw(m_prec*2) << std::left << "k"
			<< std::setw(m_prec*2) << std::left << "l"
			<< std::setw(m_prec*2) << std::left << "E"
			<< std::setw(m_prec*2) << std::left << "w"
			<< std::setw(m_prec*2) << std::left << "w_sf1"
			<< std::setw(m_prec*2) << std::left << "w_sf2"
			<< std::setw(m_prec*2) << std::left << "w_nsf"
			<< std::endl;

		for(t_size i=0; i<num_qs; ++i)
		{
			t_real h = std::lerp(h_start, h_end, t_real(i)/t_real(num_qs-1));
			t_real k = std::lerp(k_start, k_end, t_real(i)/t_real(num_qs-1));
			t_real l = std::lerp(l_start, l_end, t_real(i)/t_real(num_qs-1));

			auto energies_and_correlations = GetEnergies(h, k, l, false);
			for(const auto& E_and_S : energies_and_correlations)
			{
				ofstr
					<< std::setw(m_prec*2) << std::left << h
					<< std::setw(m_prec*2) << std::left << k
					<< std::setw(m_prec*2) << std::left << l
					<< std::setw(m_prec*2) << E_and_S.E
					<< std::setw(m_prec*2) << E_and_S.weight
					<< std::setw(m_prec*2) << E_and_S.weight_channel[0]
					<< std::setw(m_prec*2) << E_and_S.weight_channel[1]
					<< std::setw(m_prec*2) << E_and_S.weight_channel[2]
					<< std::endl;
			}
		}
	}


	/**
	 * load a configuration from a file
	 */
	bool Load(const std::string& filename)
	{
		// properties tree
		boost::property_tree::ptree node;

		// read xml file
		std::ifstream ifstr{filename};
		boost::property_tree::read_xml(ifstr, node);

		// check signature
		if(auto optInfo = node.get_optional<std::string>("magdyn.meta.info");
			!optInfo || !(*optInfo==std::string{"magdyn_tool"}))
		{
			return false;
		}

		const auto &magdyn = node.get_child("magdyn");
		return Load(magdyn);
	}


	/**
	 * save a configuration to a file
	 */
	bool Save(const std::string& filename) const
	{
		// properties tree
		boost::property_tree::ptree node;

		// write signature
		node.put<std::string>("meta.info", "magdyn_tool");
		node.put<std::string>("meta.date",
			tl2::epoch_to_str<t_real>(tl2::epoch<t_real>()));

		boost::property_tree::ptree root_node;
		root_node.put_child("magdyn", node);

		// write xml file
		std::ofstream ofstr{filename};
		if(!ofstr)
			return false;

		ofstr.precision(m_prec);
		boost::property_tree::write_xml(ofstr, node,
			boost::property_tree::xml_writer_make_settings(
				'\t', 1, std::string{"utf-8"}));
		return true;
	}


	/**
	 * load a configuration from a property tree
	 */
	bool Load(const boost::property_tree::ptree& node)
	{
		Clear();

		// variables
		if(auto vars = node.get_child_optional("variables"); vars)
		{
			for(const auto &var : *vars)
			{
				auto name = var.second.get_optional<std::string>("name");
				if(!name)
					continue;

				Variable variable;
				variable.name = *name;
				variable.value = var.second.get<t_cplx>("value", 0.);

				AddVariable(std::move(variable));
			}
		}

		// atom sites
		if(auto sites = node.get_child_optional("atom_sites"); sites)
		{
			for(const auto &site : *sites)
			{
				AtomSite atom_site;

				atom_site.name = site.second.get<std::string>("name", "n/a");
				atom_site.index = m_sites.size();

				atom_site.pos = tl2::create<t_vec_real>(
				{
					site.second.get<t_real>("position_x", 0.),
					site.second.get<t_real>("position_y", 0.),
					site.second.get<t_real>("position_z", 0.),
				});

				atom_site.spin_dir[0] = site.second.get<std::string>("spin_x", "0"),
				atom_site.spin_dir[1] = site.second.get<std::string>("spin_y", "0"),
				atom_site.spin_dir[2] = site.second.get<std::string>("spin_z", "1"),

				atom_site.spin_mag = site.second.get<t_real>("spin_magnitude", 1.);
				atom_site.g = -2. * tl2::unit<t_mat>(3);

				AddAtomSite(std::move(atom_site), false);
			}
		}

		// exchange terms
		if(auto terms = node.get_child_optional("exchange_terms"); terms)
		{
			for(const auto &term : *terms)
			{
				ExchangeTerm exchange_term;

				exchange_term.name = term.second.get<std::string>("name", "n/a");
				exchange_term.index = m_exchange_terms.size();
				exchange_term.atom1 = term.second.get<t_size>("atom_1_index", 0);
				exchange_term.atom2 = term.second.get<t_size>("atom_2_index", 0);

				exchange_term.dist = tl2::create<t_vec_real>(
				{
					term.second.get<t_real>("distance_x", 0.),
					term.second.get<t_real>("distance_y", 0.),
					term.second.get<t_real>("distance_z", 0.),
				});

				exchange_term.J = term.second.get<std::string>("interaction", "0");

				exchange_term.dmi[0] = term.second.get<std::string>("dmi_x", "0");
				exchange_term.dmi[1] = term.second.get<std::string>("dmi_y", "0");
				exchange_term.dmi[2] = term.second.get<std::string>("dmi_z", "0");

				AddExchangeTerm(std::move(exchange_term), false);
			}
		}

		// external field
		if(auto field = node.get_child_optional("field"); field)
		{
			ExternalField thefield;

			thefield.mag = 0.;
			thefield.align_spins = false;

			thefield.dir = tl2::create<t_vec_real>(
			{
				field->get<t_real>("direction_h", 0.),
				field->get<t_real>("direction_k", 0.),
				field->get<t_real>("direction_l", 1.),
			});

			if(auto optVal = field->get_optional<t_real>("magnitude"))
				thefield.mag = *optVal;

			if(auto optVal = field->get_optional<bool>("align_spins"))
				thefield.align_spins = *optVal;

			SetExternalField(thefield);
		}

		// temperature
		m_temperature = node.get<t_real>("temperature", -1.);

		// ordering vector
		if(auto ordering = node.get_child_optional("ordering"); ordering)
		{
			t_vec_real ordering_vec = tl2::create<t_vec_real>(
			{
				ordering->get<t_real>("h", 0.),
				ordering->get<t_real>("k", 0.),
				ordering->get<t_real>("l", 0.),
			});

			SetOrderingWavevector(ordering_vec);
		}

		// rotation axis
		if(auto axis = node.get_child_optional("rotation_axis"); axis)
		{
			t_vec_real rotaxis = tl2::create<t_vec_real>(
			{
				axis->get<t_real>("h", 1.),
				axis->get<t_real>("k", 0.),
				axis->get<t_real>("l", 0.),
			});

			SetRotationAxis(rotaxis);
		}

		CalcAtomSites();
		CalcExchangeTerms();
		return true;
	}


	/**
	 * save a configuration to a property tree
	 */
	bool Save(boost::property_tree::ptree& node) const
	{
		// external field
		node.put<t_real>("field.direction_h", m_field.dir[0]);
		node.put<t_real>("field.direction_k", m_field.dir[1]);
		node.put<t_real>("field.direction_l", m_field.dir[2]);
		node.put<t_real>("field.magnitude", m_field.mag);
		node.put<bool>("field.align_spins", m_field.align_spins);

		// ordering vector
		if(m_ordering.size() == 3)
		{
			node.put<t_real>("ordering.h", m_ordering[0]);
			node.put<t_real>("ordering.k", m_ordering[1]);
			node.put<t_real>("ordering.l", m_ordering[2]);
		}

		// rotation axis
		if(m_rotaxis.size() == 3)
		{
			node.put<t_real>("rotation_axis.h", m_rotaxis[0]);
			node.put<t_real>("rotation_axis.k", m_rotaxis[1]);
			node.put<t_real>("rotation_axis.l", m_rotaxis[2]);
		}

		// temperature
		node.put<t_real>("temperature", m_temperature);

		// variables
		for(const auto& var : GetVariables())
		{
			boost::property_tree::ptree itemNode;
			itemNode.put<std::string>("name", var.name);
			itemNode.put<t_cplx>("value", var.value);

			node.add_child("variables.variable", itemNode);
		}

		// atom sites
		for(const auto& site : GetAtomSites())
		{
			boost::property_tree::ptree itemNode;
			itemNode.put<std::string>("name", site.name);
			itemNode.put<t_real>("position_x", site.pos[0]);
			itemNode.put<t_real>("position_y", site.pos[1]);
			itemNode.put<t_real>("position_z", site.pos[2]);
			itemNode.put<std::string>("spin_x", site.spin_dir[0]);
			itemNode.put<std::string>("spin_y", site.spin_dir[1]);
			itemNode.put<std::string>("spin_z", site.spin_dir[2]);
			itemNode.put<t_real>("spin_magnitude", site.spin_mag);

			node.add_child("atom_sites.site", itemNode);
		}

		// exchange terms
		for(const auto& term : GetExchangeTerms())
		{
			boost::property_tree::ptree itemNode;
			itemNode.put<std::string>("name", term.name);
			itemNode.put<t_size>("atom_1_index", term.atom1);
			itemNode.put<t_size>("atom_2_index", term.atom2);
			itemNode.put<t_real>("distance_x", term.dist[0]);
			itemNode.put<t_real>("distance_y", term.dist[1]);
			itemNode.put<t_real>("distance_z", term.dist[2]);
			itemNode.put<std::string>("interaction", term.J);
			itemNode.put<std::string>("dmi_x", term.dmi[0]);
			itemNode.put<std::string>("dmi_y", term.dmi[1]);
			itemNode.put<std::string>("dmi_z", term.dmi[2]);

			node.add_child("exchange_terms.term", itemNode);
		}

		return true;
	}


protected:
	/**
	 * converts the rotation matrix rotating the local spins to ferromagnetic
	 * [001] directions into the vectors comprised of the matrix columns
	 * @see equation (9) from (Toth 2015)
	 */
	std::tuple<t_vec, t_vec> R_to_uv(const t_mat& R)
	{
		t_vec u = tl2::col<t_mat, t_vec>(R, 0)
			+ s_imag * tl2::col<t_mat, t_vec>(R, 1);
		t_vec v = tl2::col<t_mat, t_vec>(R, 2);

		return std::make_tuple(u, v);
	}


	/**
	 * spin rotation of equation (9) from (Toth 2015)
	 * rotate local spin to ferromagnetic [001] direction
	 */
	std::tuple<t_vec, t_vec> spin_to_uv(const t_vec& spin_dir)
	{
		auto [spin_re, spin_im] =
			tl2::split_cplx<t_vec, t_vec_real>(spin_dir);

		if(!tl2::equals_0<t_vec_real>(spin_im, m_eps))
			std::cerr << "Spin vector should be purely real." << std::endl;

		// only use real part, imaginary part should be zero
		//spin_re /= tl2::norm<t_vec_real>(spin_re);
		t_mat_real _rot = tl2::rotation<t_mat_real, t_vec_real>(
			spin_re, m_zdir, &m_rotaxis, m_eps);

		t_mat rot = tl2::convert<t_mat, t_mat_real>(_rot);
		return R_to_uv(rot);
	}


private:
	// atom sites
	std::vector<AtomSite> m_sites{};
	std::vector<AtomSiteCalc> m_sites_calc{};

	// magnetic couplings
	std::vector<ExchangeTerm> m_exchange_terms{};
	std::vector<ExchangeTermCalc> m_exchange_terms_calc{};

	// open variables
	std::vector<Variable> m_variables{};

	// external field
	ExternalField m_field{};
	// matrix to rotate field into the [001] direction
	t_mat m_rot_field = tl2::unit<t_mat>(3);

	// ordering wave vector for incommensurate structures
	t_vec_real m_ordering = tl2::zero<t_vec_real>(3);
	t_vec_real m_rotaxis = tl2::create<t_vec_real>({1., 0., 0.});

	// direction to rotation spins into, usually [001]
	t_vec_real m_zdir = tl2::create<t_vec_real>({0., 0., 1.});

	// temperature (-1: disable bose factor)
	t_real m_temperature{-1};

	// bose cutoff energy to avoid infinities
	t_real m_bose_cutoff{0.025};

	// settings
	bool m_is_incommensurate{false};
	bool m_force_incommensurate{false};
	bool m_unite_degenerate_energies{true};

	// settings for cholesky decomposition
	t_size m_tries_chol{50};
	t_real m_delta_chol{0.01};

	// precisions
	t_real m_eps{1e-6};
	int m_prec{6};

	// constants
	static constexpr const t_cplx s_imag {t_real(0), t_real(1)};
	static constexpr const t_real s_twopi {t_real(2)*tl2::pi<t_real>};
};

}
#endif
