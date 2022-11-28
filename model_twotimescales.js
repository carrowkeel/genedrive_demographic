
const defaults = () => ({
	'target_steps': 100,
	's': '0.6',
	'c': '0.8',
	'h': '0.5',
	'q1': '0.5',
	'q2': '0',
	'm': '0.01',
	'n1': '1',
	'n2': '1',
	'r': '1',
	'gametic': 1
});

const step = (params, _step, t) => {
	const variables = Object.assign({}, params, _step);
	variables['sn'] = 0.5 * (1 - variables['c']) * (1 - variables['h'] * variables['s']);
	variables['sc'] = variables['gametic'] === 1 ? variables['c'] * (1 - variables['h'] * variables['s']) : variables['c'] * (1 - variables['s']);
const sc = params.gametic === 1 ? params.c * (1 - params.h * params.s) : params.c * (1 - params.s);
	variables['qm1'] = !variables['n1'] && !variables['n2'] ? variables['q1'] : (!variables['n1'] ? variables['q2'] : (variables['q1'] * (1 - variables['m']) + variables['q2'] * variables['m'] * variables['n2']/variables['n1']) * 1/(1 - variables['m'] + variables['m'] * variables['n2']/variables['n1']));
	variables['qm2'] = !variables['n1'] && !variables['n2'] ? variables['q2'] : (!variables['n2'] ? variables['q1'] : (variables['q2'] * (1 - variables['m']) + variables['q1'] * variables['m'] * variables['n1']/variables['n2']) * 1/(1 - variables['m'] + variables['m'] * variables['n1']/variables['n2']));
	variables['w1'] = variables['qm1']**2 * (1 - variables['s']) + 2 * variables['qm1'] * (1 - variables['qm1']) * (2 * variables['sn'] + variables['sc']) + (1 - variables['qm1'])**2;
	variables['w2'] = variables['qm2']**2 * (1 - variables['s']) + 2 * variables['qm2'] * (1 - variables['qm2']) * (2 * variables['sn'] + variables['sc']) + (1 - variables['qm2'])**2;
	variables['q1'] = (variables['qm1']**2 * (1 - variables['s']) + 2 * variables['qm1'] * (1 - variables['qm1']) * (variables['sn'] + variables['sc'])) / variables['w1'];
	variables['q2'] = (variables['qm2']**2 * (1 - variables['s']) + 2 * variables['qm2'] * (1 - variables['qm2']) * (variables['sn'] + variables['sc'])) / variables['w2'];
	variables['n1'] = 1 - variables['q1']**variables['d'];
	variables['n2'] = 1 - variables['q2']**variables['d'];
	return {'q1': variables['q1'], 'q2': variables['q2'], 'n1': variables['n1'], 'n2': variables['n2']};
};

const sum = (arr) => arr.reduce((a,v)=>a+v,0);
const mean = (arr) => arr.length === 0 ? 0 : arr.reduce((a,v)=>a+v,0)/arr.length;
const range = (start,end) => Array.from(Array(end-start)).map((v,i)=>i+start);

const outcome = (params, steps) => {
	const freq = [steps[steps.length-1].q1, steps[steps.length-1].q2];
	const population = steps.reduce((a,step) => range(0, 2).map(i => a[i].concat(step[`n${i+1}`])), [[], []]);
	const max_values = population.map(p => p.filter((v,i) => i > 0 && i < p.length-1 && ((p[i-1] < v && v > p[i+1]) || (p[i-1] > v && v < p[i+1]))));
	const min = population.map(pop => Math.min.apply(null, pop));
	const pop_end = population.map(v => v[v.length-1]);
	const gen_threshold = 0.5;
	const pop_threshold = 0.9;
	switch(true) {
		case min[0] < pop_threshold && min[1] < pop_threshold:
			return 'collapse';
		case freq[0] >= gen_threshold && freq[1] >= gen_threshold:
			return 'spillover';
		case freq[0] < gen_threshold && freq[1] < gen_threshold && min[0] < pop_threshold && min[1] >= pop_threshold && pop_end[0] >= pop_threshold && pop_end[1] >= pop_threshold && population[0][0] > min[0]:
			return 'rescue';
		case freq[0] < gen_threshold && freq[1] < gen_threshold && pop_end[0] >= pop_threshold && pop_end[1] >= pop_threshold:
			return 'failure';
		case freq[1] < gen_threshold && min[0] < pop_threshold && min[1] >= pop_threshold && pop_end[0] < pop_threshold && pop_end[1] >= pop_threshold:
			return 'suppression';
		case freq[0] >= gen_threshold && freq[1] < gen_threshold && pop_end[0] >= pop_threshold && pop_end[1] >= pop_threshold:
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

const result = (params, steps) => {
	const model_outcome = outcome(params, steps);
	return {
		outcome: model_outcome
	};
};

const run = (params) => {
	const repeats = range(0, params.repeats || 1).map(i => sim(params));
	if (repeats.length === 1)
		return repeats[0];
	return Object.keys(repeats[0]).reduce((a,stat) => {
		const values = repeats.map(result => result.stat);
		return Object.assign(a, {[stat]: stat === 'outcome' ?
			common_value(values) :
			mean(values)});
	}, {});
};

export {defaults, step, result, run}
