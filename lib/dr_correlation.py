"""
Module that describes Duns & Ros correlation
"""
import math as mt
from typing import Optional

import lib.constants as cnst
import lib.friction as fr


class DunsRos:
    """
    Class for calculating pressure gradient in wells using Duns & Ros correlation
    """

    def __init__(self):
        self.fric = fr.Friction()
        self.const = cnst.dr_const
        self.fp = None
        self.vsl = None
        self.vsg = None
        self.vsm = None
        self.n_d = None
        self.rho_s_kgm3 = None
        self.ek = None
        self.a = None
        self.dp_dl_grav = None
        self.dp_dl_fric = None
        self.dp_dl = None

    @staticmethod
    def calc_fp(n_gv: float, n_b_s: float, n_s_tr: float, n_tr_m: float) -> int:
        """
        Definition of flow pattern
        :param n_gv: gas velocity number, dimensionless
        :param n_b_s: bubble-slug flows boundary, dimensionless
        :param n_s_tr: slug-transition flows boundary, dimensionless
        :param n_tr_m: transition-mist flows boundary, dimensionless
        :return: flow pattern, dimensionless
            * 0 - Bubble, plug and part of the froth flows pattern
            * 1 - Slug flow pattern
            * 2 - Transition flow pattern
            * 3 - Mist flow pattern
        """
        if n_gv <= n_b_s:
            fp = 0
        elif n_b_s < n_gv <= n_s_tr:
            fp = 1
        elif n_s_tr < n_gv < n_tr_m:
            fp = 2
        else:
            fp = 3
        return fp

    def calc_hl(
            self,
            rho_lrc_kgm3: float,
            sigma_l_nm: float,
            n_gv: float,
            n_lv: float,
            n_l: float,
    ) -> float:
        """
        Calculating the liquid holdup
        :param rho_lrc_kgm3: liquid's density in P,T conditions, kg/m3
        :param sigma_l_nm: surface tension factor, N/m
        :param n_gv: gas velocity number, dimensionless
        :param n_lv: liquid velocity number, dimensionless
        :param n_l: liquid viscosity number, dimensionless
        :return: liquid hold-up, dimensionless
        """
        if self.fp == 0:
            f1 = self.const["f1"](n_l)
            f2 = self.const["f2"](n_l)
            f3 = self.const["f3"](n_l)
            f4 = self.const["f4"](n_l)
            f3_hatch = f3 - f4 / self.n_d
            s = f1 + f2 * n_lv + f3_hatch * (n_gv / (1 + n_lv)) ** 2
        else:
            f5 = self.const["f5"](n_l)
            f6 = self.const["f6"](n_l)
            f7 = self.const["f7"](n_l)
            f6_hash = 0.029 * self.n_d + f6
            s = (1 + f5) * ((n_gv ** 0.982 + f6_hash) / (1 + f7 * n_lv) ** 2)
        # Slip velocity
        vs = s / (rho_lrc_kgm3 / (9.81 * sigma_l_nm)) ** 0.25
        # Liquid hold-up
        hl = (vs - self.vsm + ((self.vsm - vs) ** 2 + 4 * vs * self.vsl) ** 0.5) / (2 * vs)
        return hl

    @staticmethod
    def calc_rho_s(
            rho_lrc_kgm3: float,
            rho_grc_kgm3: float,
            fp: int,
            hl: float,
            ll: float,
            n_gv: float,
            n_tr_m: float,
            ek: float,
            a: float,
    ) -> float:
        """
        Calculation of the density of the mixture with slip effect
        :param rho_lrc_kgm3: liquid's density in P,T conditions, kg/m3
        :param rho_grc_kgm3: density of the gas in P,T conditions, kg/m3
        :param fp: flow pattern, dimensionless
            * 0 - Bubble, plug and part of the froth flows pattern
            * 1 - Slug flow pattern
            * 2 - Transition flow pattern
            * 3 - Mist flow pattern
        :param hl: liquid hold-up, dimensionless
        :param ll: liquid load-up (gas-liquid velocity ratio), dimensionless
        :param n_gv: gas velocity number, dimensionless
        :param n_tr_m: transition-mist flows boundary, dimensionless
        :param ek:  kinetic energy, dimensionless
        :param a: transition's coefficient, dimensionless
        :return: density of the mixture with slip effect, kg/m3
        """
        rho_s_1 = rho_lrc_kgm3 * hl + rho_grc_kgm3 * (1 - hl)
        if fp == 2:
            rho_grc_new = rho_grc_kgm3 * n_gv / n_tr_m
        else:
            rho_grc_new = rho_grc_kgm3
        rho_s_2 = (rho_lrc_kgm3 * ll + rho_grc_new * (1 - ll)) / (1 - ek)
        rho_s_kgm3 = rho_s_1 * a + rho_s_2 * (1 - a)
        return rho_s_kgm3

    @staticmethod
    def calc_grav(rho_s_kgm3: float, theta_deg: float, c_grav: Optional[float] = 1) -> float:
        """
        Calculation of the gravity pressure gradient
        :param rho_s_kgm3: density of the mixture with slip effect, kg/m3
        :param theta_deg: pipe angle, deg
        :param c_grav: calibration factor for gravity, dimensionless
        :return: gravity pressure gradient, Pa/m
        """
        dp_dl_grav = (rho_s_kgm3 * 9.81 * mt.sin(theta_deg / 180 * mt.pi) * c_grav)
        return dp_dl_grav

    def calc_ff_l(self, d_tub: float, roughness: float, rho_lrc_kgm3: float, mul_rc_cp: float) -> float:
        """
        Calculation of the friction factor of the liquid phase by correlation Duns & Ros (fp = 0, 1)
        :param d_tub: diameter of the tube, m
        :param roughness: tube's roughness, m
        :param rho_lrc_kgm3: liquid's density in P, T conditions, kg/m3
        :param mul_rc_cp: viscosity of the liquid in P, T conditions , sPs
        :return: friction factor of the liquid phase, dimensionless
        """
        n_re_l = self.fric.calc_n_re(rho_lrc_kgm3, self.vsl, mul_rc_cp, d_tub)
        f1 = self.fric.calc_ff(n_re_l, d_tub, roughness)
        f2 = self.const["f2_y"](f1 / 4 * self.vsg / self.vsl * self.n_d ** (2 / 3))
        f3 = 1 + f1 / 4 * (self.vsg / (50 * self.vsl)) ** 0.5
        return f1 * f2 / f3

    def calc_ff_g(
        self,
        d_tub: float,
        roughness: float,
        rho_lrc_kgm3: float,
        rho_grc_kgm3: float,
        mul_rc_cp: float,
        mug_rc_cp: float,
        sigma_l_nm: float,
        vsg: float
    ) -> float:
        """
        Calculation of the friction factor of the gas phase by correlation Duns & Ros (fp = 3)
        :param d_tub: diameter of the tube, m
        :param roughness: tube's roughness, m
        :param rho_lrc_kgm3: liquid's density in P, T conditions, kg/m3
        :param rho_grc_kgm3: gas density in P, T conditions, kg/m3
        :param mul_rc_cp: viscosity of the liquid in P, T conditions , sPs
        :param mug_rc_cp: viscosity of the gas in P, T conditions , sPs
        :param sigma_l_nm: surface tension factor, N/m
        :param vsg: gas velocity, m/s
        :return: friction factor of the gas phase, dimensionless
        """
        # Weber's number
        n_we = rho_grc_kgm3 * vsg ** 2 * roughness / sigma_l_nm
        # Dimensionless interaction between viscosity and surface tension
        n_mu = (mul_rc_cp * 1e-3) ** 2 / rho_lrc_kgm3 / sigma_l_nm / roughness
        if n_mu * n_we <= 0.005:
            rel_rough_new = 0.0749 * sigma_l_nm / rho_grc_kgm3 / vsg ** 2 / d_tub
        else:
            rel_rough_new = (
                0.3713
                * sigma_l_nm
                * (n_mu * n_we) ** 0.302
                / rho_grc_kgm3
                / vsg ** 2
                / d_tub
            )
        if rel_rough_new > 0.05:
            ff_g = 4 * (
                (4 * mt.log10(0.27 * rel_rough_new)) ** (-2)
                + 0.067 * rel_rough_new ** 1.73
            )
        else:
            n_re_g = self.fric.calc_n_re(rho_grc_kgm3, vsg, mug_rc_cp, d_tub)
            ff_g = self.fric.calc_ff(n_re_g, d_tub, roughness)
        return ff_g

    def calc_fric(
            self,
            d_tub: float,
            roughness: float,
            rho_lrc_kgm3: float,
            rho_grc_kgm3: float,
            mul_rc_cp: float,
            mug_rc_cp: float,
            sigma_l_nm: float,
            c_fric: Optional[float] = 1
    ) -> float:
        """
        Calculation of the friction pressure gradient
        :param d_tub: diameter of the tube, m
        :param roughness: tube's roughness, m
        :param rho_lrc_kgm3: liquid's density in P,T conditions, kg/m3
        :param rho_grc_kgm3: gas density in P,T conditions, kg/m3
        :param mul_rc_cp: viscosity of the liquid in P, T conditions , sPs
        :param mug_rc_cp: viscosity of the gas in P, T conditions , sPs
        :param sigma_l_nm: surface tension factor, N/m
        :param c_fric: calibration factor for friction, dimensionless
        :return: friction pressure gradient, Pa/m
        """
        ff_l = self.calc_ff_l(d_tub, roughness, rho_lrc_kgm3, mul_rc_cp)
        grad_fric_1 = ff_l * rho_lrc_kgm3 * self.vsl * self.vsm / (2 * d_tub)
        if self.fp == 2 or self.fp == 3:
            d_new = d_tub - 2 * roughness
            vsg_new = self.vsg * d_tub ** 2 / d_new ** 2
            ff_g = self.calc_ff_g(
                d_new,
                roughness,
                rho_lrc_kgm3,
                rho_grc_kgm3,
                mul_rc_cp,
                mug_rc_cp,
                sigma_l_nm,
                vsg_new
            )
            grad_fric_2 = ff_g * rho_grc_kgm3 * vsg_new ** 2 / (2 * d_new) / (1 - self.ek)
        else:
            grad_fric_2 = 0
        dp_dl_fr = (grad_fric_1 * self.a + grad_fric_2 * (1 - self.a)) * c_fric
        return dp_dl_fr
    
    def calc_params(
            self,
            d_tub: float,
            ql_rc_m3s: float,
            qg_rc_m3s: float,
            rho_lrc_kgm3: float,
            rho_grc_kgm3: float,
            mul_rc_cp: float,
            sigma_l_nm: float,
            p: float
    ) -> None:
        """
        Calculating preparatory parameters of Duns & Ros correlation
        :param d_tub: diameter of the tube, m
        :param ql_rc_m3s: rate of the liquid, m3/s
        :param qg_rc_m3s: rate of the gas, m3/s
        :param rho_lrc_kgm3: liquid's density in P,T conditions, kg/m3
        :param rho_grc_kgm3: gas density in P,T conditions, kg/m3
        :param mul_rc_cp: viscosity of the liquid in P, T conditions , sPs
        :param sigma_l_nm: surface tension factor, N/m
        :param p: pressure on the current stage, Pa
        :return: none
        """
        # Velocities of the liquid, gas and mixture
        self.vsl = ql_rc_m3s / (mt.pi * d_tub ** 2 / 4)
        self.vsg = qg_rc_m3s / (mt.pi * d_tub ** 2 / 4)
        self.vsm = self.vsl + self.vsg
        # Dimensionless parameters
        n_gv = self.vsg * (rho_lrc_kgm3 / (9.81 * sigma_l_nm)) ** 0.25
        n_lv = self.vsl * (rho_lrc_kgm3 / (9.81 * sigma_l_nm)) ** 0.25
        self.n_d = d_tub * (rho_lrc_kgm3 * 9.81 / sigma_l_nm) ** 0.5
        n_l = (mul_rc_cp * 1e-3) * (9.81 / (rho_lrc_kgm3 * sigma_l_nm ** 3)) ** 0.25
        # Factors for the bubble-slug flows boundary
        l1 = self.const["l1"](self.n_d)
        l2 = self.const["l2"](self.n_d)
        # Boundaries between flow regimes
        n_b_s = l1 + l2 * n_lv
        n_s_tr = 50 + 36 * n_lv
        n_tr_m = 75 + 84 * n_lv ** 0.75
        # Flow pattern
        self.fp = self.calc_fp(n_gv, n_b_s, n_s_tr, n_tr_m)
        # Interpolation factor for the transition regime (fp = 2)
        if self.fp == 0 or self.fp == 1:
            self.a = 1
        elif self.fp == 2:
            self.a = (n_tr_m - n_gv) / (n_tr_m - n_s_tr)
        else:
            self.a = 0
        # Liquid load-up (gas-liquid velocity ratio)
        ll = max(ql_rc_m3s / (ql_rc_m3s + qg_rc_m3s), 0.000001)
        # Mixtures density by gas-liquid velocity ratio
        rho_n_kgm3 = rho_lrc_kgm3 * ll + rho_grc_kgm3 * (1 - ll)
        # Kinetic energy
        if self.fp == 2 or self.fp == 3:
            self.ek = self.vsm * self.vsg * rho_n_kgm3 / p if p != 0 else 0
        else:
            self.ek = 0
        # Liquid hold-up
        if self.fp != 3:
            hl = self.calc_hl(rho_lrc_kgm3, sigma_l_nm, n_gv, n_lv, n_l)
            if self.fp != 2:
                hl = hl
            else:
                hl = hl * self.a + ll * (1 - self.a)
        else:
            hl = ll
        # Density of the mixture with slip effect
        self.rho_s_kgm3 = self.calc_rho_s(rho_lrc_kgm3, rho_grc_kgm3, self.fp, hl, ll, n_gv, n_tr_m, self.ek, self.a)

    def calc_grad(
            self,
            d_tub: float,
            roughness: float,
            theta_deg: float,
            ql_rc_m3s: float,
            qg_rc_m3s: float,
            rho_lrc_kgm3: float,
            rho_grc_kgm3: float,
            mul_rc_cp: float,
            mug_rc_cp: float,
            sigma_l_nm: float,
            p: float
    ) -> float:
        """
        Calculation of the pressure gradient
        :param d_tub: diameter of the tube, m
        :param roughness: tube's roughness, m
        :param theta_deg: pipe angle, deg
        :param ql_rc_m3s: rate of the liquid, m3/s
        :param qg_rc_m3s: rate of the gas, m3/s
        :param rho_lrc_kgm3: liquid's density in P,T conditions, kg/m3
        :param rho_grc_kgm3: gas density in P,T conditions, kg/m3
        :param mul_rc_cp: viscosity of the liquid in P, T conditions , sPs
        :param mug_rc_cp: viscosity of the gas in P, T conditions , sPs
        :param sigma_l_nm: surface tension factor, N/m
        :param p: pressure on the current stage, Pa
        :return: pressure gradient, Pa/m
        """
        # Preparatory parameters
        self.calc_params(d_tub, ql_rc_m3s, qg_rc_m3s, rho_lrc_kgm3, rho_grc_kgm3, mul_rc_cp, sigma_l_nm, p)
        # Gravity pressure gradient
        self.dp_dl_grav = self.calc_grav(self.rho_s_kgm3, theta_deg)
        # Friction pressure gradient
        self.dp_dl_fric = self.calc_fric(d_tub, roughness, rho_lrc_kgm3, rho_grc_kgm3, mul_rc_cp, mug_rc_cp, sigma_l_nm)
        # Total pressure gradient
        self.dp_dl = self.dp_dl_grav + self.dp_dl_fric
        return self.dp_dl
