from .skin import Skin
import math as m

class PerfDW:
    """
    Класс для расчета перфорированной наклонно-направленной скважины

    full_perf_s - метод расчета полностью перфорированной скважины

    part_perf_s - метод расчета частично перфорированной скважины
    """
    def __init__(self) -> None:
        self.skin = Skin()

    def full_perf_s(self, k: float, kd: float, rw: float, rd: float, teta: float,
                    phi: float, Lp: float, rp: float, ns: float, kh: float, kv: float,
                    kcz: float, rcz: float, model: int, h: float, hw: float, zw: float) -> float:
        """
        Метод расчета полностью перфорированной скважины

        Parameters
        ----------
        :param k: начальная проницаемость, мД;
        :param kd: измененная проницаемость, мД;
        :param rw: радиус скважины, м;
        :param rd: радиус зоны с проницаемостью, измененной по сравнению с проницаемостью пласта, м;
        :param teta: угол отклонения ствола скважины от вертикали, градусы
        :param phi: фазировка перфорационных зарядов, градусы;
        :param Lp: длина перфорационных каналов, м;
        :param rp: радиус перфорационных каналов, м;
        :param ns: плотность перфорационных отверстий, отв/м;
        :param kh: проницаемость пласта в латеральном направлении, мД;
        :param kv: проницаемость пласта в вертикальном направлении, мД;
        :param kcz: проницаемость зоны разрушения породы вокруг перфорационных каналов, мД;
        :param rcz: радиус зоны разрушения породы вокруг перфорационных каналов, м;
        :param model: 0 - Корреляция Cinco-Ley; 1 - Корреляция Ozkan-Raghavan;
        :param h: мощность пласта, м;
        :param hw: мощность вскрытого интервала, открытого для притока (0<=hw<=h), м;
        :param zw: расстояние от подошвы пласта до центра интервала, открытого для притока (hw/2<=zw<=h-hw/2), м;

        ----------
        """
        St = m.cos(teta)*(self.skin.calc_Sd(k, kd, rw, rd) + k/kd*self.skin.calc_Sp(phi, rw, Lp, rp, ns, kh, kv) + self.skin.calc_Scz(ns, Lp, k, kcz, kd, rcz, rp)) + self.skin.calc_Steta(model, teta, kh, kv, h, hw, rw, zw)
        return St


    def part_perf_s(self, k: float, kd: float, rw: float, rd: float, teta: float,
                    phi: float, Lp: float, rp: float, ns: float, kh: float, kv: float,
                    kcz: float, rcz: float, h: float, hw: float, zw: float) -> float:
        """
        Метод расчета частично перфорированной скважины

        Parameters
        ----------
        :param k: начальная проницаемость, мД;
        :param kd: измененная проницаемость, мД;
        :param rw: радиус скважины, м;
        :param rd: радиус зоны с проницаемостью, измененной по сравнению с проницаемостью пласта, м;
        :param teta: угол отклонения ствола скважины от вертикали, градусы
        :param phi: фазировка перфорационных зарядов, градусы;
        :param Lp: длина перфорационных каналов, м;
        :param rp: радиус перфорационных каналов, м;
        :param ns: плотность перфорационных отверстий, отв/м;
        :param kh: проницаемость пласта в латеральном направлении, мД;
        :param kv: проницаемость пласта в вертикальном направлении, мД;
        :param kcz: проницаемость зоны разрушения породы вокруг перфорационных каналов, мД;
        :param rcz: радиус зоны разрушения породы вокруг перфорационных каналов, м;
        :param h: мощность пласта, м;
        :param hw: мощность вскрытого интервала, открытого для притока (0<=hw<=h), м;
        :param zw: расстояние от подошвы пласта до центра интервала, открытого для притока (hw/2<=zw<=h-hw/2), м;

        ----------
        """
        St = m.cos(teta)*(self.skin.calc_Sd(k, kd, rw, rd) + k/kd*self.skin.calc_Sp(phi, rw, Lp, rp, ns, kh, kv) + self.skin.calc_Scz(ns, Lp, k, kcz, kd, rcz, rp)) + self.skin.calc_Sopp(teta, kh, kv, h, hw, rw, zw)
        return St
