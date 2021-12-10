from django.shortcuts import render
import numpy as np
import json
from ninja import NinjaAPI
from app.skin.skin import *
from app.skin.uncased_vertical_well import UncasedVW
from app.skin.perforated_vertical_well import PerfVW
from app.skin.unanchored_directional_well import UnanchDW
from app.skin.perforated_directional_well import PerfDW

api = NinjaAPI()

def index(request):
	"""
	# Метод представления главной страницы

	"""
	return render(request, "app/index.html")

@api.post("/plot0")
def plot0(request, data0):
	data0 = dict_verify(json.loads(data0))
	uncased_vertical_well = UncasedVW()
	perf_vertical_well = PerfVW()
	unanchored_directional_well = UnanchDW()
	perforated_directional_well = PerfDW()
	if data0['type'] == 10:
		S = uncased_vertical_well.perfect_s(data0['k'], data0['kd'], data0['rw'], data0['rd'])
		res_l = 'Производительность совершенной по степени вскрытия вертикальной скважины'
	elif data0['type'] == 11:
		S = uncased_vertical_well.unperfect_s(data0['model'], data0['h'], data0['hw'], data0['rw'], data0['zw'], data0['k'], data0['kv'], data0['k'], data0['kd'], data0['rd'], data0['y'])
		res_l = 'Производительность несовершенной по степени вскрытия вертикальной скважины'
	elif data0['type'] == 20:
		S = perf_vertical_well.full_perf_s(data0['k'], data0['kd'], data0['rw'], data0['rd'], data0['phi'], data0['Lp'], data0['rp'], data0['ns'], data0['k'], data0['kv'], data0['kcz'], data0['rcz'])
		res_l = 'Производительность полностью перфорированной вертикальной скважины'
	elif data0['type'] == 21:
		S = perf_vertical_well.part_perf_s(data0['k'], data0['kd'], data0['rw'], data0['rd'], data0['phi'], data0['Lp'], data0['rp'], data0['ns'], data0['k'], data0['kv'], data0['kcz'], data0['rcz'], data0['h'], data0['hw'], data0['model'], data0['zw'], data0['y'])
		res_l = 'Производительность частично перфорированной вертикальной скважины'
	elif data0['type'] == 30:
		S = unanchored_directional_well.perfect_s(data0['k'], data0['kd'], data0['rw'], data0['rd'], data0['model'], data0['teta'], data0['k'], data0['kv'], data0['h'], data0['hw'], data0['zw'])
		res_l = 'Производительность совершенной по степени вскрытия необсаженной наклонно-направленной скважины'
	elif data0['type'] == 31:
		S = unanchored_directional_well.unperfect_s(data0['k'], data0['kd'], data0['rw'], data0['rd'], data0['h'], data0['Lwpc'], data0['teta'], data0['k'], data0['kv'], data0['hw'], data0['zw'])
		res_l = 'Производительность несовершенной по степени вскрытия наклонно-направленной скважины'
	elif data0['type'] == 40:
		S = perforated_directional_well.full_perf_s(data0['k'], data0['kd'], data0['rw'], data0['rd'], data0['teta'], data0['phi'], data0['Lp'], data0['rp'], data0['ns'], data0['k'], data0['kv'], data0['kcz'], data0['rcz'], data0['model'], data0['h'], data0['hw'], data0['zw'])
		res_l = 'Производительность полностью перфорированной наклонно-направленной скважины'
	elif data0['type'] == 41:
		S = perforated_directional_well.part_perf_s(data0['k'], data0['kd'], data0['rw'], data0['rd'], data0['teta'], data0['phi'], data0['Lp'], data0['rp'], data0['ns'], data0['k'], data0['kv'], data0['kcz'], data0['rcz'], data0['h'], data0['hw'], data0['zw'])
		res_l = 'Производительность частично перфорированной наклонно-направленной скважины'

	
	q = q_well(data0['k'], data0['h'], data0['Pres'], data0['Pwf'], data0['mu'], data0['B'], data0['re'], data0['rw'], S)
	r_arr, p_arr = np.linspace(data0['rw'], 100, 500).tolist(), []
	for r in r_arr:
		p_arr.append(p_ss_atma(
			data0['Pres'],
			q,
			data0['mu'],
			data0['B'],
			data0['k'],
			data0['h'],
			data0['re'],
			S,
			r
			))
	return {"res": f'{res_l} - {round(q,1)} м3/сут', "r_arr":r_arr , "p_arr":p_arr }

def dict_verify(dict_: dict):
	"""
	## Функция валидации необязательных параметров
	"""
	for key in dict_:
		if isinstance(dict_[key], str):
			if dict_[key] != '':
				dict_[key] = float(dict_[key])
			else:
				dict_[key] = None
	return dict_