"""
Module that describes single-phase flow friction by Moody
"""
import math as mt


class Friction:
    """
    Class for calculating the single-phase flow friction factor
    """

    @staticmethod
    def calc_n_re(
            rho_n: float,
            v_ms: float,
            mu: float,
            d_tub: float
    ) -> float:
        """
        Function for calculating the Reynold's number
        :param rho_n: density of the current phase, kg/m3
        :param v_ms: velocity of the current phase, m3/s
        :param mu: viscosity of the current phase, sPs
        :param d_tub: diameter of the tube, m
        :return: Reynold's number, dimensionless
        """
        return rho_n * v_ms * d_tub / mu * 1e3

    @staticmethod
    def calc_ff(
            n_re: float,
            d_tub: float,
            roughness: float,
    ) -> float:
        """
        Function for calculating the coefficient of friction by Churchill correlation
        :param n_re: Reynold's number, dimensionless
        :param d_tub: diameter of the tube, m
        :param roughness: tube's roughness, m
        :return: friction factor, dimensionless
        """
        a = (-2.457 * mt.log((7 / n_re) ** 0.9 + 0.27 * (roughness / d_tub))) ** 16
        b = (37530 / n_re) ** 16
        ff = 8 * ((8 / n_re) ** 12 + 1 / (a + b) ** 1.5) ** (1 / 12)
        return ff
