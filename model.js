const randint = (m,m1,g) => Math.floor(g() * (m1 - m)) + m;
const round = (n,p) => { var f = Math.pow(10, p); return Math.round(n * f) / f };
const sum = (arr) => arr.reduce((a,v)=>a+v,0);
const mean = (arr) => arr.length === 0 ? 0 : arr.reduce((a,v)=>a+v,0)/arr.length;
const range = (start,end) => Array.from(Array(end-start)).map((v,i)=>i+start);
const cumsum = arr => { let sum = 0; const out = []; for (const v of arr) { sum += v; out.push(sum) } return out; }
const common_value = values => Object.entries(values.reduce((a,value) => Object.assign(a, {[value]: (a[value] || 0) + 1}), {})).sort((a,b)=>b[1]-a[1])[0][0];

const initial = params => {
	return {q_1: params.q_1, q_2: params.q_2, p_1: 0, p_2: 0, N_1: params.K - params.K * params.eps, N_2: params.K - params.K * params.eps};
};

const step = (params, _step = initial(params), t) => {
	if (t === 0)
		return _step;
	const {q_1, q_2, p_1, p_2, N_1, N_2} = _step;
	const sn = 0.5 * (1 - params.c) * (1 - params.h * params.s);
	const sc = params.gametic === 1 ? params.c * (1 - params.h * params.s) : params.c * (1 - params.s);
	const qm = [
		((1 - params.m) * q_1 + params.m * (N_2/N_1) * q_2) * 1/(1 - params.m + params.m * (N_2/N_1)),
		((1 - params.m) * q_2 + params.m * (N_1/N_2) * q_1) * 1/(1 - params.m + params.m * (N_1/N_2))
	];
	const pm = [
		((1 - params.m) * p_1 + params.m * (N_2/N_1) * p_2) * 1/(1 - params.m + params.m * (N_2/N_1)),
		((1 - params.m) * p_2 + params.m * (N_1/N_2) * p_1) * 1/(1 - params.m + params.m * (N_1/N_2))
	];
	const n = [(1 - params.m) * N_1 + params.m * N_2, (1 - params.m) * N_2 + params.m * N_1];
	if (n[0] < 0 || n[1] < 0)
		return false;
	// Selection
	const w = qm.map((q,i) => 2*q*(1-q-pm[i])*(2*sn + sc)+2*pm[i]*(1-q-pm[i])*(1-params.rc/2)+2*pm[i]*q*(1-params.s*params.h-params.rc/2)+q**2*(1-params.s)+(1-q-pm[i])**2+pm[i]**2*(1-params.rc));
	const qn = qm.map((q,i) => (q*(1-q-pm[i])*2*(sn+sc)+q*pm[i]*(1-params.s*params.h-params.rc/2)+q**2*(1-params.s))/w[i]);
	const pn = qm.map((q,i) => (pm[i]*(1-q-pm[i])*(1-params.rc/2)+pm[i]*q*(1-params.s*params.h-params.rc/2)+pm[i]**2*(1-params.rc))/w[i]);
	// Mutation
	qn.forEach((q,i) => {
		if (pn[i] === 0 && params.r_freq > 0 && Math.random() < params.r_freq*2*_n[i]*((1-q)**2 + q*(1-q)*(1-params.c)))
			pn[i] += 0.001;
		if (pn[i] === 0 && params.nhej_freq > 0 && Math.random() < params.nhej_freq*2*_n[i]*params.c*q*(1-q))
			pn[i] += 0.001;
	});
	const n_tag = n.map((N,i) => {
		const r0 = params.r0;
		const r1 = r0 - (1 - Math.E**(-params.s * params.d)) * (params.s + r0);
		const rh = r0 * (1 - params.h) + params.h * r1;
		const r = (1-qn[i])**2 * r0 + 2*qn[i]*(1-qn[i])*((1-params.c)*rh + params.c*(params.gametic ? rh : r1)) + qn[i]**2 * r1;
		return N + r * N * (1 - (N / params.K));
	});
	return (isNaN(qn[0]) || isNaN(qn[1])) ? false : {q_1: qn[0], q_2: qn[1], p_1: pn[0], p_2: pn[1], N_1: n_tag[0], N_2: n_tag[1]};
};

const outcome = (params, steps) => {
	const freq = [steps[steps.length-1].q_1, steps[steps.length-1].q_2];
	const resistance_freq = [steps[steps.length-1].p_1, steps[steps.length-1].p_2];
	const population = steps.reduce((a,step) => range(0, 2).map(i => a[i].concat(step[`N_${i+1}`] / params.K)), [[], []]);
	const max_values = population.map(p => p.filter((v,i) => i > 0 && i < p.length-1 && ((p[i-1] < v && v > p[i+1]) || (p[i-1] > v && v < p[i+1]))));
	const non_converging = population.map(p => p.filter((v,i) => i > 0 && i < p.length-1 && p[i-1] < v && v > p[i+1]).filter((v,i,arr)=>(i>0?v-arr[i-1]:-1)>=0).length > 0);
	const oscillating = max_values.map(p => mean(p.map((v,i,arr)=>i>0?Math.abs(v-arr[i-1]):0).slice(1)) > 0.25);
	const min = population.map(pop => Math.min.apply(null, pop));
	const pop_end = population.map(v => v[v.length-1]);
	const gen_threshold = 0.5;
	const pop_threshold = [0.9, 0.5, 0.2];
	switch(true) {
		case min[0] < pop_threshold[0] && min[1] < pop_threshold[0]: // freq[0] >= gen_threshold && freq[1] >= gen_threshold && 
			return 'collapse';
		case freq[0] >= gen_threshold && freq[1] >= gen_threshold:
			return 'spillover';
		case freq[0] < gen_threshold && freq[1] < gen_threshold && (resistance_freq[0] >= 0.1 || resistance_freq[1] >= 0.1) && min[0] >= pop_threshold[0] && min[1] >= pop_threshold[0] && pop_end[0] >= pop_threshold[0] && pop_end[1] >= pop_threshold[0]:
			return 'res_failure';
		case freq[0] < gen_threshold && freq[1] < gen_threshold && (resistance_freq[0] >= 0.1 || resistance_freq[1] >= 0.1) && min[0] < pop_threshold[0] && min[1] >= pop_threshold[0] && pop_end[0] >= pop_threshold[0] && pop_end[1] >= pop_threshold[0]:
			return 'res_rescue';
		case freq[0] < gen_threshold && freq[1] < gen_threshold && min[0] < pop_threshold[0] && min[1] >= pop_threshold[0] && pop_end[0] >= pop_threshold[0] && pop_end[1] >= pop_threshold[0]:
			return 'rescue';
		case freq[0] < gen_threshold && freq[1] < gen_threshold && min[0] >= pop_threshold[0] && min[1] >= pop_threshold[0] && pop_end[0] >= pop_threshold[0] && pop_end[1] >= pop_threshold[0]:
			return 'failure';
		case freq[0] > 0 && freq[1] < gen_threshold && (oscillating[0] || non_converging[0]):
			return 'oscillating';
		case freq[1] < gen_threshold && min[0] < pop_threshold[0] && min[1] >= pop_threshold[0] && pop_end[0] < pop_threshold[0] && pop_end[1] >= pop_threshold[0]:
			return 'suppression';
		case freq[0] >= gen_threshold && freq[1] < gen_threshold && pop_end[0] >= pop_threshold[0] && pop_end[1] >= pop_threshold[0]:
			return 'dte';
		default:
			return 'undefined';
	}
};

const sim = (params) => {
	const storage = [];
	let s = undefined;
	while (s !== false && storage.length <= params.target_steps) {
		s = step(params, s, storage.length);
		if (s)
			storage.push(s);
	}
	return result(params, storage);
};

const mcrit = (params) => {
	let m = 0;
	let incr = -2;
	const precision = -6;
	while (true) {
		const result = sim(Object.assign({}, params, {m, stat: 'outcome'}));
		if (!['failure', 'spillover', 'collapse'].includes(result.outcome)) {
			m += Math.pow(10, incr);
		} else if (m !== 0) {
			m -= Math.pow(10, incr);
			if (incr === precision)
				return {mcrit: Math.min(1, m)};
			incr--;
			m += Math.pow(10, incr);
		} else {
			return {mcrit: Math.min(1, m)};
		}
	}
};

const defaults = () => ({
	s: 0.5,
	c: 0.8,
	h: 0.5,
	d: 5,
	r0: 1,
	m: 0.01,
	q_1: 0.8,
	q_2: 0,
	r_freq: 0,
	nhej_freq: 0,
	rc: 0,
	K: 1,
	gametic: 0,
	eps: 0.01,
	target_steps: 100,
	repeats: 1,
	outcome: 'suppression'
});

const run = (params) => {
	const repeats = range(0, params.repeats).map(i => (params.stat === 'mcrit' ? mcrit : sim)(params));
	if (repeats.length === 1)
		return repeats[0];
	return Object.keys(repeats[0]).reduce((a,stat) => {
		const values = repeats.map(result => result[stat]);
		return Object.assign(a, {[stat]: stat === 'outcome' ?
			values.filter(outcome => outcome === params.outcome).length / values.length :
			mean(values)});
	}, {});
};

const result = (params, steps) => {
	const model_outcome = outcome(params, steps);
	const non_target = steps.sort((a,b) => b.q_2 - a.q_2)[0].q_2;
	const population = steps.reduce((a,step) => range(0, 2).map(i => a[i].concat(step[`N_${i+1}`] / params.K)), [[], []]);
	return {
		outcome: model_outcome,
		nontarget: non_target,
		suppression: ['failure', 'rescue', 'oscillating', 'suppression'].includes(model_outcome) ? 1 - Math.min.apply(null, population[0]) : 0,
		suppression_sum: ['failure', 'rescue', 'oscillating', 'suppression'].includes(model_outcome) 
 ? sum(population[0].map(v => 1 - v)) / population[0].length : 0,
		suppression_duration: ['failure', 'rescue', 'oscillating', 'suppression'].includes(model_outcome) ? Math.min(1, population[0].filter(v => v<0.9).length / population[0].length) : 0
	};
};

export { step, run, defaults }
