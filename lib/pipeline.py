"""
Module for calculation of pressure and temperature in the pipe
"""
import numpy as np
from numpy.typing import ArrayLike
from scipy.integrate import solve_ivp

import lib.dr_correlation as dr_cor


class Pipeline:
    """
    Class of pipe for calculating VLP curve
    """

    def __init__(self):
        self.dr = dr_cor.DunsRos()

    def integr_func(
            self,
            h: float,
            pt: tuple,
            d_tub: float,
            roughness: float,
            theta_deg: float,
            ql_rc_m3s: float,
            qg_rc_m3s: float,
            rho_lrc_kgm3: float,
            rho_grc_kgm3: float,
            mul_rc_kgm3: float,
            mug_rc_kgm3: float,
            sigma_l_nm: float,
            temp_grad: float
    ) -> tuple:
        """
        Function for pipe integration
        :param h: current depth, m
        :param pt: current pressure and temperature, Pa, K
        :param d_tub: diameter of the tube, m
        :param roughness: tube's roughness, m
        :param theta_deg: pipe angle, deg
        :param ql_rc_m3s: rate of the liquid, m3/s
        :param qg_rc_m3s: rate of the gas, m3/s
        :param rho_lrc_kgm3: liquid's density in P,T conditions, kg/m3
        :param rho_grc_kgm3: gas density in P,T conditions, kg/m3
        :param mul_rc_kgm3: viscosity of the liquid in P, T conditions , sPs
        :param mug_rc_kgm3: viscosity of the gas in P, T conditions , sPs
        :param sigma_l_nm: surface tension factor, N/m
        :param temp_grad: temperature gradient, K/100m
        :return: pressure and temperature gradients, Pa/m, K/m
        """
        p, t = pt
        dp_dl = self.dr.calc_grad(
            d_tub,
            roughness,
            theta_deg,
            ql_rc_m3s,
            qg_rc_m3s,
            rho_lrc_kgm3,
            rho_grc_kgm3,
            mul_rc_kgm3,
            mug_rc_kgm3,
            sigma_l_nm,
            p
        )
        dt_dl = temp_grad / 100
        return dp_dl, dt_dl

    def calc_pipe(
            self,
            p_wh: float,
            t_wh: float,
            h0: float,
            md_vdp: tuple,
            d_tub: float,
            roughness: float,
            theta_deg: float,
            q_fluid_m3s: float,
            gor: float,
            rho_lrc_kgm3: float,
            rho_grc_kgm3: float,
            mul_rc_kgm3: float,
            mug_rc_kgm3: float,
            sigma_l_nm: float,
            temp_grad: float
    ) -> tuple:
        """
        Calculating pressure and temperature inside the pipe
        :param p_wh: wellhead pressure, Pa
        :param t_wh: wellhead temperature, K
        :param h0: current depth, m
        :param md_vdp: measured depth, m
        :param d_tub: diameter of the tube, m
        :param roughness: tube's roughness, m
        :param theta_deg: pipe angle, deg
        :param q_fluid_m3s: rate of the fluid, m3/s
        :param gor: gas oil ratio, m3/m3
        :param rho_lrc_kgm3: liquid's density in P,T conditions, kg/m3
        :param rho_grc_kgm3: gas density in P,T conditions, kg/m3
        :param mul_rc_kgm3: viscosity of the liquid in P, T conditions , sPs
        :param mug_rc_kgm3: viscosity of the gas in P, T conditions , sPs
        :param sigma_l_nm: surface tension factor, N/m
        :param temp_grad: temperature gradient, K/100m
        :return: pressure and temperature, Pa, K
        """
        qg_rc_m3s = gor * q_fluid_m3s
        pipe = solve_ivp(
            self.integr_func,
            t_span=(h0, md_vdp),
            y0=[p_wh, t_wh],
            method='RK23',
            args=(
                d_tub,
                roughness,
                theta_deg,
                q_fluid_m3s,
                qg_rc_m3s,
                rho_lrc_kgm3,
                rho_grc_kgm3,
                mul_rc_kgm3,
                mug_rc_kgm3,
                sigma_l_nm,
                temp_grad
            )
        )
        return pipe.y[0, :], pipe.y[1, :]

    def calc_pwf(
            self,
            p_wh: float,
            t_wh: float,
            h0: float,
            md_vdp: tuple,
            d_tub: float,
            roughness: float,
            theta_deg: float,
            q_fluid_m3s: float,
            gor: float,
            rho_lrc_kgm3: float,
            rho_grc_kgm3: float,
            mul_rc_kgm3: float,
            mug_rc_kgm3: float,
            sigma_l_nm: float,
            temp_grad: float
    ) -> float:
        """
        Calculating well floor pressure
        :param p_wh: wellhead pressure, Pa
        :param t_wh: wellhead temperature, K
        :param h0: starting depth, m
        :param md_vdp: measured depth, m
        :param d_tub: diameter of the tube, m
        :param roughness: tube's roughness, m
        :param theta_deg: pipe angle, deg
        :param q_fluid_m3s: rate of the fluid, m3/s
        :param gor: gas oil ratio, m3/m3
        :param rho_lrc_kgm3: liquid's density in P,T conditions, kg/m3
        :param rho_grc_kgm3: gas density in P,T conditions, kg/m3
        :param mul_rc_kgm3: viscosity of the liquid in P, T conditions , sPs
        :param mug_rc_kgm3: viscosity of the gas in P, T conditions , sPs
        :param sigma_l_nm: surface tension factor, N/m
        :param temp_grad: temperature gradient, K/100m
        :return: well floor pressure, Pa
        """
        return self.calc_pipe(
            p_wh,
            t_wh,
            h0,
            md_vdp,
            d_tub,
            roughness,
            theta_deg,
            q_fluid_m3s,
            gor,
            rho_lrc_kgm3,
            rho_grc_kgm3,
            mul_rc_kgm3,
            mug_rc_kgm3,
            sigma_l_nm,
            temp_grad
        )[0][-1]

    def calc_vlp(
            self,
            p_wh: float,
            t_wh: float,
            h0: float,
            md_vdp: float,
            d_tub: float,
            roughness: float,
            theta_deg: float,
            q_fluid_m3s: tuple,
            gor: float,
            rho_lrc_kgm3: float,
            rho_grc_kgm3: float,
            mul_rc_kgm3: float,
            mug_rc_kgm3: float,
            sigma_l_nm: float,
            temp_grad: float
    ) -> tuple[ArrayLike, ArrayLike]:
        """
        Calculating VLP
        :param p_wh: wellhead pressure, Pa
        :param t_wh: wellhead temperature, K
        :param h0: starting depth, m
        :param md_vdp: measured depth, m
        :param d_tub: diameter of the tube, m
        :param roughness: tube's roughness, m
        :param theta_deg: pipe angle, deg
        :param q_fluid_m3s: rate of the fluid, m3/s
        :param gor: gas oil ratio, m3/m3
        :param rho_lrc_kgm3: liquid's density in P,T conditions, kg/m3
        :param rho_grc_kgm3: gas density in P,T conditions, kg/m3
        :param mul_rc_kgm3: viscosity of the liquid in P, T conditions , sPs
        :param mug_rc_kgm3: viscosity of the gas in P, T conditions , sPs
        :param sigma_l_nm: surface tension factor, N/m
        :param temp_grad: temperature gradient, K/100m
        :return:
        """
        p_wf_func = np.vectorize(self.calc_pwf)
        p_wfs = p_wf_func(
            p_wh,
            t_wh,
            h0,
            md_vdp,
            d_tub,
            roughness,
            theta_deg,
            q_fluid_m3s,
            gor,
            rho_lrc_kgm3,
            rho_grc_kgm3,
            mul_rc_kgm3,
            mug_rc_kgm3,
            sigma_l_nm,
            temp_grad
        ) / 101325
        return q_fluid_m3s, p_wfs
