from .skin import Skin


class UncasedVW:
    """
    Класс для расчета необсаженной вертикальной скважины

    calc_perfect_s - метод расчета скин-фактора совершенной скважины по степени вскрытия

    calc_unperfect_s - метод расчета скин-фактора несовершенной скважины по степени вскрытия
    """
    def __init__(self) -> None:
        self.skin = Skin()

    def perfect_s(self, k: float, kd: float, rw: float, rd: float) -> float:
        """
        Метод расчета скин-фактора совершенной скважины по степени вскрытия

        Parameters
        ----------
        :param k: начальная проницаемость, мД;
        :param kd: измененная проницаемость, мД;
        :param rw: радиус скважины, м;
        :param rd: радиус  зоны  с  проницаемостью,  измененной  по  сравнению  с проницаемостью пласта, м;

        ----------
        """
        St = self.skin.calc_Sd(k, kd, rw, rd)
        return St

    def unperfect_s(self, model: int, h: float, hw: float, rw: float, zw: float, kh: float, kv: float, k: float, kd: float, rd: float, y: float = 1) -> float:
        """
        Метод расчета скин-фактора несовершенной скважины по степени вскрытия

        Parameters
        ----------
        :param model: 0 - Корреляция Papatzacos; 1 - Корреляция Vrbik;
        :param h: мощность пласта, м;
        :param hw: мощность вскрытого интервала, открытого для притока (0<=hw<=h), м;
        :param rw: радиус скважины, м;
        :param zw: расстояние от подошвы пласта до центра интервала, открытого для притока (hw/2<=zw<=h-hw/2), м;
        :param kh: проницаемость пласта в латеральном направлении, мД;
        :param kv: проницаемость пласта в вертикальном направлении, мД;
        :param k: начальная проницаемость, мД;
        :param kd: измененная проницаемость, мД;
        :param rd: радиус  зоны  с  проницаемостью,  измененной  по  сравнению  с проницаемостью пласта, м;
        :param y: коэффициент несовершенства степени вскрытия;

        ----------
        """
        St = 1/y*h/hw*self.skin.calc_Sd(k, kd, rw, rd) + self.skin.calc_Spp(model, h, hw, rw,zw, kh, kv)
        return St
