import math as m
import numpy as np

def p_ss_atma(p_res_atma = 250,
              q_liq_sm3day = 50,
              mu_cP = 1,
              B_m3m3 = 1.2,
              k_mD = 40,
              h_m = 10,
              r_e = 240,
              S = 0,
              r = 0.1):
    """
    функция расчета давления в произвольной точке пласта для стационарного решения 
    уравнения фильтрации 
    p_res_atma - пластовое давление, давление на контуре питания
    q_liq_sm3day - дебит жидкости на поверхности в стандартных условиях
    mu_cP - вязкость нефти (в пластовых условиях)
    B_m3m3 - объемный коэффициент нефти 
    k_mD - проницаемость пласта
    h_m - мощность пласта
    r_e - радиус контрура питания 
    S - скин фактора (расчетный)
    r - расстояние на котором проводится расчет
    """
    return p_res_atma - 18.41 * q_liq_sm3day*mu_cP*B_m3m3/k_mD/h_m * (np.log(r_e/r)+S-0.75)

def q_well(k: float, h: float, Pres: float, Pwf: float, mu: float, B: float, re: float, rw: float, S: float):
    """
    Функция для расчета производительности скважины на псевдоустановившемся режиме

    Parameters
    ----------
    :param k: проницаемость пласта, мД
    :param h: толщина пласта, м
    :param Pres: среднее пластовое давление, атм
    :param Pwf: забойное давление, атм
    :param mu: вякость нефти, сПз
    :param B: объемный коэффициент, м3/м3
    :param re: радиус контура, м
    :param rw: радиус скважины, м
    :param S: скин фактора (расчетный)

    :return q: производительность скважины, м3/сут

    ----------
    """
    return k*h*(Pres-Pwf)/(18.4*mu*B*(m.log1p(re/rw)+S-0.75))
class Skin:
    def calc_Sd(self, k: float, kd: float, rw: float, rd: float) -> float:
        """
        Метод расчета механического скин-фактора
        
        Parameters
        ----------
        :param k: начальная проницаемость, мД;
        :param kd: измененная проницаемость, мД;
        :param rw: радиус скважины, м;
        :param rd: радиус  зоны  с  проницаемостью,  измененной  по  сравнению  с проницаемостью пласта, м;

        ----------
        """
        return (k/kd - 1)*m.log1p(rd/rw)

    def calc_Spp(self, model: int, h: float, hw: float, rw: float, zw: float, kh: float, kv: float) -> float:
        """
        Метод расчета скин-фактора за счет частичного вскрытия
        
        Parameters
        ----------
        :param model: 0 - Корреляция Papatzacos; 1 - Корреляция Vrbik;
        :param h: мощность пласта, м;
        :param hw: мощность вскрытого интервала, открытого для притока (0<=hw<=h), м;
        :param rw: радиус скважины, м;
        :param zw: расстояние от подошвы пласта до центра интервала, открытого для притока (hw/2<=zw<=h-hw/2), м;
        :param kh: проницаемость пласта в латеральном направлении, мД;
        :param kv: проницаемость пласта в вертикальном направлении, мД;

        ----------
        """
        hwh = hw/h
        hd = h/rw*(kh/kv)**0.5
        if model == 0:
            ls = ((hwh)/(2+hwh))*(((zw+hw/4)*(h/zw+hw/4))/((zw-hw/4)*(h-zw-hw/4)))**0.5
            Spp = (1/hwh-1)*m.log1p(3.14*hd/2) + m.log1p(ls)/hwh
            return Spp
        else:
            f_0 = self.Vrbik_func(0, hd)
            f_1 = self.Vrbik_func(hwh, hd)
            f_2 = self.Vrbik_func(2*zw/h, hd)
            f_3 = self.Vrbik_func((2*zw+hw)/h, hd)
            f_4 = self.Vrbik_func((2*zw-hw)/h, hd)
            Spp = (1/hwh-1)*(1.2704+m.log1p(hd))-(1/hwh)**2*(f_0-f_1+f_2-0.5*f_3-0.5*f_4)
            return Spp

    def Vrbik_func(self, y: float, hd: float) -> float:
        if y == 2 or y == 0:
            return 2*m.log1p(2) + 1/(3.14*hd)*m.log1p(0.1053/hd**2)
        else:
            return y*m.log1p(y) + (2-y)*m.log1p(2-y)+1/(3.14*hd)*m.log1p(m.sin(3.14*y/2)**2+0.1053/hd**2)

    def calc_Sp(self, phi: float, rw: float, Lp: float,
        rp: float, ns: float, kh: float, kv: float) -> float:
        """
        Метод расчета скин-фактора за счет перфорации по корреляции Karakas-Tariq
        
        Parameters
        ----------
        :param phi: фазировка перфорационных зарядов, градусы;
        :param rw: радиус скважины, м;
        :param Lp: длина перфорационных каналов, м;
        :param rp: радиус перфорационных каналов, м;
        :param ns: плотность перфорационных отверстий, отв/м;
        :param kh: проницаемость пласта в латеральном направлении, мД;
        :param kv: проницаемость пласта в вертикальном направлении, мД;

        ----------
        """
        self.set_coeff(phi)
        Sp = self.calc_Sh(rw, self.a, Lp) + self.calc_Sv(
            [self.a1, self.a2, self.b1, self.b2],
            rp, ns, Lp, kh, kv
        ) + self.calc_Swb(
            [self.c1, self.c2],
            rw, Lp
        )
        return Sp

    def set_coeff(self, phi: int) -> None:
        """
        Метод определения числовых коэффициентов, зависящих от фазировки перфорационных зарядов

        Parameters
        ----------
        :param phi: фазировка перфорационных зарядов, градусы;

        ----------
        """
        self.table = {
            0: [0.250, -2.091, 0.0453, 5.1313, 1.8672, 1.6/10**1, 2.675],
            360: [0.250, -2.091, 0.0453, 5.1313, 1.8672, 1.6/10**1, 2.675],
            180: [0.500, -2.025, 0.0943, 3.0373, 1.8115, 2.6/10**2, 4.532],
            120: [0.648, -2.018, 0.0634, 1.6136, 1.7770, 6.6/10**3, 5.320],
            90: [0.726, -1.905, 0.1038, 1.5674, 1.6935, 1.9/10**3, 6.155],
            60: [0.813, -1.898, 0.1023, 1.3654, 1.6490, 3.0/10**4, 7.509],
            45: [0.860, -1.788, 0.2398, 1.1915, 1.6392, 4.6/10**5, 8.791],
        }
        self.a = self.table[phi][0]
        self.a1 = self.table[phi][1]
        self.a2 = self.table[phi][2]
        self.b1 = self.table[phi][3]
        self.b2 = self.table[phi][4]
        self.c1 = self.table[phi][5]
        self.c2 = self.table[phi][6]

    def calc_Sh(self, rw: float, a: float, Lp: float) -> float:
        """
        Метод расчета скин-фактора за счет схождения потока к перфорационным каналам в горизонатльной плоскости
        
        Parameters
        ----------
        :param rw: радиус скважины, м;
        :param a: набор численных коэффициентов, зависящих от фазировки перфорационных каналов;
        :param Lp: длина перфорационных каналов, м;

        ----------
        """
        # rwe - эффективный радиус скважины с учетом длины перфорационных каналов, м
        if a == 0:
            rwe = Lp/4
        else:
            rwe = a*(rw+Lp)
        Sh = m.log1p(rw/rwe)
        return Sh

    def calc_Sv(self, coef_list: list, rp: float, ns: float, Lp: float, kh: float, kv: float) -> float:
        """
        Метод расчета скин-фактора за счет схождения потока к перфорационным каналам в вертикальной плоскости

        Parameters
        ----------
        :param coef_list = [a1, a2, b1, b2]: набор числовых констант, зависящие от фазировки перфорационных зарядов;
        :param rp: радиус перфорационных каналов, м;
        :param ns: плотность перфорационных отверстий, отв/м;
        :param Lp: длина перфорационных каналов, м;
        :param kh: проницаемость пласта в латеральном направлении, мД;
        :param kv: проницаемость пласта в вертикальном направлении, мД;

        ----------
        """
        dzp = 1/ns # расстояние между перфорационными отверстиями, м
        rpd = rp/(2*dzp)*(1+(kv/kh)**0.5)
        a = coef_list[0]*m.log1p(rpd) + coef_list[1]
        b = coef_list[2]*rpd + coef_list[3]
        zpd = dzp/Lp*(kh/kv)**0.5
        Sv = 10**a*zpd**(b-1)*rpd**b
        return Sv

    def calc_Swb(self, coef_list: list, rw: float, Lp: float) -> float:
        """
        Метод расчета скин-фактора за счет самого ствола скважины

        Parameters
        ----------
        :param coef_list = [с1, с2]: набор числовых констант, зависящие от фазировки перфорационных зарядов;
        :param rw: радиус скважины, м;
        :param Lp: длина перфорационных каналов, м;

        ----------
        """
        rwd = rw/(rw+Lp)
        Swb = coef_list[0]*m.exp(coef_list[1]*rwd)
        return Swb

    def calc_Scz(self, ns: float, Lp: float, k: float, kcz: float, kd: float, rcz: float, rp: float) -> float:
        """
        Метод расчет скин-фактора за счет зоны разрушения овркуг перфорационных каналов

        Parameters
        ----------
        :param ns: плотность перфорационных отверстий, отв/м;
        :param Lp: длина перфорационных каналов, м;
        :param k: начальная проницаемость, мД;
        :param kcz: проницаемость зоны разрушения породы вокруг перфорационных каналов, мД;
        :param kd: измененная проницаемость, мД;
        :param rcz: радиус зоны разрушения породы вокруг перфорационных каналов, м;
        :param rp: радиус перфорационных каналов, м;

        ----------
        """
        dzp = 1/ns # расстояние между перфорационными отверстиями, м
        Scz = dzp/Lp*(k/kcz-k/kd)*m.log1p(rcz/rp)
        return Scz

    def calc_Steta(self, model: int, teta: float, kh: float, kv: float, h: float, hw: float, rw: float, zw: float) -> float:
        """
        Метод расчета геометрического скин-фактора за счет отклонения скважины от вертикали, определяемый 
        по корреляциям Cinco-Ley или Ozkan-Raghavan для скважины, полностью вскрывающей продуктивный пласт

        Parameters
        ----------
        :param model: 0 - Корреляция Cinco-Ley; 1 - Корреляция Ozkan-Raghavan;
        :param teta: угол отклонения ствола скважины отв вертикали, градусы;
        :param kh: проницаемость пласта в латеральном направлении, мД;
        :param kv: проницаемость пласта в вертикальном направлении, мД;
        :param h: мощность пласта, м;
        :param hw: мощность вскрытого интервала, открытого для притока (0<=hw<=h), м;
        :param rw: радиус скважины, м;
        :param zw: расстояние от подошвы пласта до центра интервала, открытого для притока (hw/2<=zw<=h-hw/2), м;

        ----------
        """
        if model == 0:
            teta_ = m.atan((kv/kh)*m.tan(teta))
            hd = hw/rw*(kh/kv)**0.5
            Steta = -(teta_/41)**2.06 - (teta_/56)**1.865*m.log1p(hd/100)
        else:
            Steta = self.calc_Sopp(teta, kh, kv, h, hw, rw, zw)
        return Steta

    def g_func(self, x: float, y: float, a: float, b: float) -> float:
        return 0.25*((x-b)*m.log1p((x-b)**2 + y**2) - (x-a)*m.log1p((x-a)**2 + y**2) - y/2*(m.atan((x-a)/y) - m.atan((x-b)/y)))
        
    def calc_Sopp(self, teta: float, kh: float, kv: float, h: float, hw: float, rw: float, zw: float) -> float:
        """
        Метод расчета геометрического скин-фактора за счет отклонения скважины от вертикали и за счет частичного вскрытия
        по корреляции Ozkan-Raghavan

        Parameters
        ----------
        :param teta: угол отклонения ствола скважины отв вертикали, градусы;
        :param kh: проницаемость пласта в латеральном направлении, мД;
        :param kv: проницаемость пласта в вертикальном направлении, мД;
        :param h: мощность пласта, м;
        :param hw: мощность вскрытого интервала, открытого для притока (0<=hw<=h), м;
        :param rw: радиус скважины, м;
        :param zw: расстояние от подошвы пласта до центра интервала, открытого для притока (hw/2<=zw<=h-hw/2), м;

        ----------
        """
        teta_ = m.atan((kv/kh)*m.tan(teta))
        hd = hw/rw*(kh/kv)**0.5
        hwd = hw/rw*(kh/kv*m.cos(teta)**2 + m.sin(teta)**2)**0.5
        zwd = zw/rw*(kh/kv)**0.5
        rd = (1 + 0.09*hwd**2*m.sin(teta_))**0.5
        y = m.acos((0.3*hwd*m.sin(teta_)**2)/rd)
        if zw >= h/2:
            zd = zwd + 0.3*hwd*m.cos(teta_)
        else:
            zd = zwd - 0.3*hwd*m.cos(teta_)
        e = (zd-zwd)*m.cos(teta_)**2
        yi = (3.14*rd*m.sin(y))/(hd*m.sin(teta_))
        F = -hd/(2*hwd)*(m.log1p(1 - 2*m.exp(-yi)*m.cos(3.14)*((zd + zwd + e)/hd) + m.exp(-2*yi)) + 
            m.log1p(1 - 2*m.exp(-yi)*m.cos(3.14)*((zd - zwd - e)/hd) + m.exp(-2*yi)))
        Sopp = 1 + 2/(hwd*m.sin(teta))*self.g_func(rd*m.cos(y), rd*m.sin(y), -hwd/2*m.sin(teta), hwd/2*m.sin(teta)) + F
        return Sopp