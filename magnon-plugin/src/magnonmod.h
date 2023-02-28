/**
 * S(q,w) module for magnon dynamics
 * @author Tobias Weber <tweber@ill.fr>
 * @date jan-2022
 * @license GPLv2, see 'LICENSE' file
 */

#ifndef __MAGNON_SQW_MOD_H__
#define __MAGNON_SQW_MOD_H__

#include "tools/monteconvo/sqwbase.h"
#include "tlibs/math/linalg.h"
#include "tlibs2/libs/magdyn.h"


class MagnonMod : public SqwBase
{
	public:
		using SqwBase::t_var;

		using t_size = std::size_t;
		using t_real = t_real_reso;
		using t_cplx = std::complex<t_real>;
		using t_vec = tl::ublas::vector<t_real>;

		using t_vec_real = tl2::vec<t_real, std::vector>;
		using t_mat_real = tl2::mat<t_real, std::vector>;

		using t_vec_cplx = tl2::vec<t_cplx, std::vector>;
		using t_mat_cplx = tl2::mat<t_cplx, std::vector>;

		using t_magdyn = tl2_mag::MagDyn<
			t_mat_cplx, t_vec_cplx,
			t_mat_real, t_vec_real,
			t_cplx, t_real, t_size>;


	protected:
		t_magdyn m_dyn{};

		// peak width
		t_real m_sigma = t_real(0.025);

		// S(q,E) scaling factor
		t_real m_S0 = t_real(1.);

		// incoherent amplitude and width
		t_real m_incoh_amp = t_real(0.);
		t_real m_incoh_sigma = t_real(0.025);

		// polarisation channel, -1: unpolarised
		int m_channel{-1};


	public:
		MagnonMod();
		MagnonMod(const std::string& strCfgFile);
		virtual ~MagnonMod();

		virtual std::tuple<std::vector<t_real>, std::vector<t_real>>
			disp(t_real dh, t_real dk, t_real dl) const override;
		virtual t_real operator()(t_real dh, t_real dk, t_real dl, t_real dE) const override;

		virtual std::vector<t_var> GetVars() const override;
		virtual void SetVars(const std::vector<t_var>&) override;
		virtual bool SetVarIfAvail(const std::string& strKey, const std::string& strNewVal) override;

		virtual SqwBase* shallow_copy() const override;
};

#endif
