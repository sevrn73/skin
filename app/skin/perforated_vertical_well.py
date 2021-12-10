from .skin import Skin


class PerfVW:
    """
    Класс для расчета перфорированной вертикальной скважины

    full_perf_s - метод расчета полностью перфорированной скважины

    part_perf_s - метод расчета частично перфорированной скважины
    """
    def __init__(self) -> None:
        self.skin = Skin()

    def full_perf_s(self, k: float, kd: float, rw: float, rd: float, phi: float, Lp: float,
        rp: float, ns: float, kh: float, kv: float, kcz: float, rcz: float) -> float:
        """
        Метод расчета полностью перфорированной скважины

        Parameters
        ----------
        :param k: начальная проницаемость, мД;
        :param kd: измененная проницаемость, мД;
        :param rw: радиус скважины, м;
        :param rd: радиус  зоны  с  проницаемостью,  измененной  по  сравнению  с проницаемостью пласта, м;
        :param phi: фазировка перфорационных зарядов, градусы;
        :param Lp: длина перфорационных каналов, м;
        :param rp: радиус перфорационных каналов, м;
        :param ns: плотность перфорационных отверстий, отв/м;
        :param kh: проницаемость пласта в латеральном направлении, мД;
        :param kv: проницаемость пласта в вертикальном направлении, мД;
        :param kcz: проницаемость зоны разрушения породы вокруг перфорационных каналов, мД;
        :param rcz: радиус зоны разрушения породы вокруг перфорационных каналов, м;

        ----------
        """
        St = self.skin.calc_Sd(k, kd, rw, rd) + k/kd*self.skin.calc_Sp(phi, rw, Lp, rp, ns, kh, kv) + self.skin.calc_Scz(ns, Lp, k, kcz, kd, rcz, rp)
        return St


    def part_perf_s(self, k: float, kd: float, rw: float, rd: float, phi: float, Lp: float,
        rp: float, ns: float, kh: float, kv: float, kcz: float, rcz: float, h: float, hw: float, model: int, zw: float, y: float = 1) -> float:
        """
        Метод расчета частично перфорированной скважины

        Parameters
        ----------
        :param k: начальная проницаемость, мД;
        :param kd: измененная проницаемость, мД;
        :param rw: радиус скважины, м;
        :param rd: радиус  зоны  с  проницаемостью,  измененной  по  сравнению  с проницаемостью пласта, м;
        :param phi: фазировка перфорационных зарядов, градусы;
        :param Lp: длина перфорационных каналов, м;
        :param rp: радиус перфорационных каналов, м;
        :param ns: плотность перфорационных отверстий, отв/м;
        :param kh: проницаемость пласта в латеральном направлении, мД;
        :param kv: проницаемость пласта в вертикальном направлении, мД;
        :param kcz: проницаемость зоны разрушения породы вокруг перфорационных каналов, мД;
        :param rcz: радиус зоны разрушения породы вокруг перфорационных каналов, м;
        :param h: мощность пласта, м;
        :param hw: мощность интервала перфорации, м;
        :param model: 0 - Корреляция Papatzacos; 1 - Корреляция Vrbik;
        :param zw: расстояние от подошвы пласта до центра интервала, открытого для притока (hw/2<=zw<=h-hw/2), м;
        :param y: коэффициент несовершенства степени вскрытия;

        ----------
        """
        St = h/hw/y*(self.skin.calc_Sd(k, kd, rw, rd) + k/kd*self.skin.calc_Sp(phi, rw, Lp, rp, ns, kh, kv) + self.skin.calc_Scz(ns, Lp, k, kcz, kd, rcz, rp)) + self.skin.calc_Spp(model, h, hw, rw, zw, kh, kv)
        return St