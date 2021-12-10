from .skin import Skin
import math as m


class UnanchDW:
    """
    Класс для расчета необсаженной наклонно-направленной скважины

    perfect_s - метод расчета скин-фактора совершенной скважины по степени вскрытия

    unperfect_s - метод расчета скин-фактора несовершенной скважины по степени вскрытия
    """
    def __init__(self) -> None:
        self.skin = Skin()

    def perfect_s(self, k: float, kd: float, rw: float, rd: float, model: int, teta: float,
                kh: float, kv: float, h: float, hw: float, zw: float) -> float:
        """
        Метод расчета скин-фактора совершенной скважины по степени вскрытия

        Parameters
        ----------
        :param k: начальная проницаемость, мД;
        :param kd: измененная проницаемость, мД;
        :param rw: радиус скважины, м;
        :param rd: радиус зоны с проницаемостью, измененной по сравнению с проницаемостью пласта, м;
        :param model: 0 - Корреляция Cinco-Ley; 1 - Корреляция Ozkan-Raghavan;
        :param teta: угол отклонения ствола скважины отв вертикали, градусы;
        :param kh: проницаемость пласта в латеральном направлении, мД;
        :param kv: проницаемость пласта в вертикальном направлении, мД;
        :param h: мощность пласта, м;
        :param hw: мощность вскрытого интервала, открытого для притока (0<=hw<=h), м;
        :param zw: расстояние от подошвы пласта до центра интервала, открытого для притока (hw/2<=zw<=h-hw/2), м;

        ----------
        """
        St = m.cos(teta)*self.skin.calc_Sd(k, kd, rw, rd) + self.skin.calc_Steta(model, teta, kh, kv, h, hw, rw, zw)
        return St

    def unperfect_s(self, k: float, kd: float, rw: float, rd: float, h: float, Lwpc: float,
                    teta: float, kh: float, kv: float, hw: float, zw: float) -> float:
        """
        Метод расчета скин-фактора несовершенной скважины по степени вскрытия

        Parameters
        ----------
        :param k: начальная проницаемость, мД;
        :param kd: измененная проницаемость, мД;
        :param rw: радиус скважины, м;
        :param rd: радиус зоны с проницаемостью, измененной по сравнению с проницаемостью пласта, м;
        :param h: мощность пласта, м;
        :param Lwpc: длина наклонно-направленной скважины в пределах продуктивного пласта, открытая для притока флюида, м;
        :param model: 0 - Корреляция Cinco-Ley; 1 - Корреляция Ozkan-Raghavan;
        :param teta: угол отклонения ствола скважины отв вертикали, градусы;
        :param kh: проницаемость пласта в латеральном направлении, мД;
        :param kv: проницаемость пласта в вертикальном направлении, мД;
        :param hw: мощность вскрытого интервала, открытого для притока (0<=hw<=h), м;
        :param zw: расстояние от подошвы пласта до центра интервала, открытого для притока (hw/2<=zw<=h-hw/2), м;

        ----------
        """
        St = h/Lwpc*self.skin.calc_Sd(k, kd, rw, rd) + self.skin.calc_Sopp(teta, kh, kv, h, hw, rw, zw)
        return St