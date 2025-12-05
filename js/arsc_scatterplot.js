// Path to the TSV file (relative)
const TSV_PATH = '../data/arsc_gtdb_r226rep.tsv';

// DOM elements: use the searchable inputs as the main controls
const selects = {
	domain: document.getElementById('domainSearch'),
	phylum: document.getElementById('phylumSearch'),
	class: document.getElementById('classSearch'),
	order: document.getElementById('orderSearch'),
	family: document.getElementById('familySearch'),
	genus: document.getElementById('genusSearch')
};
const ySelect = document.getElementById('ySelect');
const resetBtn = document.getElementById('resetBtn');
const downloadBtn = document.getElementById('downloadBtn');
const markerSizeInput = document.getElementById('markerSize');
const markerSizeVal = document.getElementById('markerSizeVal');

const markerAlphaInput = document.getElementById('markerAlpha');
const markerAlphaVal = document.getElementById('markerAlphaVal');

// wire up marker alpha display and input
if (markerAlphaInput && markerAlphaVal) {
	markerAlphaVal.textContent = parseFloat(markerAlphaInput.value).toFixed(2);
	markerAlphaInput.addEventListener('input', (e) => {
		markerAlphaVal.textContent = parseFloat(e.target.value).toFixed(2);
		try { update(); } catch (err) { /* update may not be defined yet */ }
	});
}

// wire up marker size display and input
if (markerSizeInput && markerSizeVal) {
	markerSizeVal.textContent = markerSizeInput.value;
	markerSizeInput.addEventListener('input', (e) => {
		markerSizeVal.textContent = e.target.value;
		// redraw with new size
		try { update(); } catch (err) { /* update may not be defined yet during initial load */ }
	});
}

let rows = []; // parsed data
let header = []; // TSV header

// Utility: parse TSV to objects
function parseTSV(text) {
	const lines = text.split('\n').filter(l => l.trim() !== '');
	const header = lines[0].split('\t');
	const data = [];
	for (let i = 1; i < lines.length; i++) {
		const cols = lines[i].split('\t');
		if (cols.length !== header.length) continue;
		const obj = {};
		for (let j = 0; j < header.length; j++) {
			obj[header[j]] = cols[j];
		}
		// parse numeric fields used for plotting
		obj['AvgResMW'] = parseFloat(obj['AvgResMW']);
		obj['GC(%)'] = parseFloat(obj['GC(%)']);
		obj['sum_len'] = parseFloat(obj['sum_len']);
		obj['N_ARSC'] = parseFloat(obj['N_ARSC']);
		obj['C_ARSC'] = parseFloat(obj['C_ARSC']);
		obj['S_ARSC'] = parseFloat(obj['S_ARSC']);
		data.push(obj);
	}
	return { header, data };
}

// Populate select options with unique sorted values
function populateSelects(data) {
	const fields = ['domain','phylum','class','order','family','genus'];
	fields.forEach(f => {
		const set = new Set();
		data.forEach(r => { if (r[f]) set.add(r[f]); });
		const arr = Array.from(set).sort();
				// populate datalist for the searchable input
				setSelectOptions(selects[f], arr);
	});
}

// Get current filters
function getFilters() {
	const f = {};
	for (const k in selects) {
		const v = selects[k].value;
		if (v) f[k] = v;
	}
	return f;
}

			// Apply filters to rows
			function filterData(rows, filters) {
				const yField = ySelect && ySelect.value ? ySelect.value : 'N_ARSC';
				const xField = 'GC(%)';
				return rows.filter(r => {
					for (const k in filters) {
						if (!r[k] || r[k] !== filters[k]) return false;
					}
					// require numeric values to plot for selected axes
					if (isNaN(r[xField]) || isNaN(r[yField])) return false;
					return true;
				});
			}

			// Draw Plotly scatter. group by phylum for categorical colors.
			function drawPlot(filtered) {
				const yField = ySelect && ySelect.value ? ySelect.value : 'N_ARSC';
				const xField = 'GC(%)';
				// decide grouping (color) level:
				// choose the deepest (most specific) selected level, then group by one level below it
				// e.g. if phylum is selected -> group by class; if none selected -> default to domain
				let groupingLevel = 'domain';
				// find deepest selected level index
				let deepest = -1;
				for (let i = 0; i < LEVELS.length; i++) {
					const lv = LEVELS[i];
					const v = selects[lv] && selects[lv].value;
					if (v) deepest = i;
				}
				if (deepest === -1) {
					groupingLevel = 'domain';
				} else if (deepest < LEVELS.length - 1) {
					groupingLevel = LEVELS[deepest + 1];
				} else {
					groupingLevel = LEVELS[deepest];
				}
				const groups = {};
				// default opacity (alpha) when no filters are applied; fallback if slider missing
				const DEFAULT_OPACITY = 1.0;
				filtered.forEach(r => {
					const key = r[groupingLevel] || 'unknown';
					if (!groups[key]) groups[key] = { x: [], y: [], text: [], size: [] };
					groups[key].x.push(r[xField]);
					groups[key].y.push(r[yField]);
					// build taxonomy info (phylum -> genus)
					const taxLevels = ['phylum','class','order','family','genus'];
					const taxHtml = taxLevels.map(t => (r[t] ? `${t}: ${r[t]}` : null)).filter(Boolean).join('<br>');
					const idLine = r['id'] ? `id: ${r['id']}` : '';
					const valLine = `${xField}: ${r[xField]}, ${yField}: ${r[yField]}`;
					groups[key].text.push(`${idLine}<br>${taxHtml}<br>${valLine}`);
					// marker size: use size from slider control (fallback to 8)
					const s = (markerSizeInput && markerSizeInput.value) ? parseInt(markerSizeInput.value) : 8;
					groups[key].size.push(s);
				});

				const traces = Object.keys(groups).map(k => ({
					x: groups[k].x,
					y: groups[k].y,
					mode: 'markers',
					type: 'scattergl',
					name: k,
					text: groups[k].text,
					hovertemplate: '%{text}<extra></extra>',
						// opacity is controlled solely by the alpha slider (or DEFAULT_OPACITY fallback)
						marker: { size: groups[k].size, sizemode: 'area', opacity: ((markerAlphaInput && markerAlphaInput.value) ? parseFloat(markerAlphaInput.value) : DEFAULT_OPACITY) }
				}));

						const layout = {
						title: yField + ' vs ' + xField + ' (color: ' + groupingLevel + ')',
						xaxis: { title: xField },
						yaxis: { title: yField },
						hovermode: 'closest',
						// place legend horizontally centered below the plot area
						legend: { orientation: 'h', x: 0.5, xanchor: 'center', y: -0.18 },
						// increase bottom margin so the legend has space
						margin: { t: 80, b: 110 }
					};

				Plotly.newPlot('plot', traces, layout, {responsive: true});
			}

// Update routine: get filters, filter rows, draw
function update() {
	const filters = getFilters();
	const filtered = filterData(rows, filters);
	drawPlot(filtered);
}

			// Hierarchy levels (top -> bottom)
			const LEVELS = ['domain','phylum','class','order','family','genus'];

			// assign data-level attributes for easier handling
			LEVELS.forEach(l => { selects[l].dataset.level = l; });

			// Helper: get unique values for a given level optionally filtered by ancestor selections
			function availableOptionsForLevel(level) {
				const idx = LEVELS.indexOf(level);
				// build ancestor filters
				const ancestorFilters = {};
				for (let i = 0; i < idx; i++) {
					const a = LEVELS[i];
					const v = selects[a].value;
					if (v) ancestorFilters[a] = v;
				}
				const set = new Set();
				rows.forEach(r => {
					let ok = true;
					for (const a in ancestorFilters) {
						if (r[a] !== ancestorFilters[a]) { ok = false; break; }
					}
					if (ok && r[level]) set.add(r[level]);
				});
				return Array.from(set).sort();
			}

			// Replace options in a select (keeps the first 'All' option)
			function setSelectOptions(sel, values) {
				// remember current selection
				const cur = sel.value;
				// populate corresponding datalist if present (inputs use datalist)
				try {
					// if sel id ends with 'Search' replace with 'List', else try other patterns
					const listId = sel.id.endsWith('Search') ? sel.id.replace('Search', 'List') : sel.id.replace('Select', 'List');
					const dl = document.getElementById(listId);
					if (dl) {
						while (dl.firstChild) dl.removeChild(dl.firstChild);
						values.forEach(v => { const o = document.createElement('option'); o.value = v; dl.appendChild(o); });
					}
				} catch (e) { /* ignore */ }
				// restore if still available
				if (cur && values.includes(cur)) sel.value = cur;
				else sel.value = '';
			}

			// Refresh all select options based on current ancestor selections
			function refreshAllOptions() {
				LEVELS.forEach((level, idx) => {
					const vals = availableOptionsForLevel(level);
					setSelectOptions(selects[level], vals);
					// also update the search input's placeholder suggestions (the datalist is filled in setSelectOptions)
					const s = selects[level];
					if (s) {
						// if current input value no longer in options, clear it to avoid mismatch
						if (s.value && !vals.includes(s.value)) s.value = '';
					}
				});
			}

			// When a lower-level select is changed, set its parents automatically when possible
					function setParentsFromChild(chLevel) {
						// Start from the changed child level and walk upwards setting parents
						const chIdx = LEVELS.indexOf(chLevel);
						let curLevel = chLevel;
						let curVal = selects[chLevel].value;
						if (!curVal) return;
						// Walk up the hierarchy: for each parent level, find a row where curLevel==curVal
						// and take its parent value. Then continue using that parent as the new curLevel.
						for (let i = chIdx - 1; i >= 0; i--) {
							const parentLevel = LEVELS[i];
							// find a row where the current level's value matches curVal and a parent exists
							const row = rows.find(r => r[curLevel] === curVal && r[parentLevel]);
							if (row) {
								selects[parentLevel].value = row[parentLevel];
								// move up: parent becomes the current level for next iteration
								curLevel = parentLevel;
								curVal = row[parentLevel];
							} else {
								// cannot determine parent for this curVal -> stop
								break;
							}
						}
					}

			// General handler when any select changes
			function onSelectChange(e) {
				const level = e.target.dataset.level;
				const idx = LEVELS.indexOf(level);
				const val = e.target.value;
				if (!val) {
					// if cleared, just refresh options below this level
					// clear deeper selections
					for (let i = idx + 1; i < LEVELS.length; i++) selects[LEVELS[i]].value = '';
					refreshAllOptions();
					update();
					return;
				}
				// If a child (lower) was changed, set parents
				// Determine if change originated at a deeper level than some currently selected parents
				// We'll always attempt to set parents from this level upwards
				setParentsFromChild(level);
				// After setting parents, refresh options to reflect consistent children sets
				refreshAllOptions();
				update();
			}

					// Wire selects (inputs)
					Object.values(selects).forEach(s => s.addEventListener('change', onSelectChange));
					if (ySelect) ySelect.addEventListener('change', update);

					// Wire search inputs (datalists) so selecting/typing a suggestion updates the control
					Object.keys(selects).forEach(level => {
						const input = selects[level];
						if (!input) return;
						input.addEventListener('change', (ev) => {
							const v = ev.target.value;
							const allowed = availableOptionsForLevel(level);
							if (v && allowed.includes(v)) {
								// set value (already set) and trigger change handling
								onSelectChange({ target: input });
							} else if (!v) {
								input.value = '';
								onSelectChange({ target: input });
							} else {
								// typed value not in allowed list -> ignore or clear
								// keep it cleared to avoid inconsistent filtering
								input.value = '';
								onSelectChange({ target: input });
							}
						});
					});
			resetBtn.addEventListener('click', () => {
				// clear taxonomy filters
				Object.values(selects).forEach(s => s.value = '');
				// reset datalist options and redraw
				refreshAllOptions();
				// reset marker size and alpha sliders to defaults if present
				if (markerSizeInput) { markerSizeInput.value = 8; }
				if (markerSizeVal) { markerSizeVal.textContent = markerSizeInput ? markerSizeInput.value : '8'; }
				if (markerAlphaInput) { markerAlphaInput.value = 1.0; }
				if (markerAlphaVal) { markerAlphaVal.textContent = markerAlphaInput ? parseFloat(markerAlphaInput.value).toFixed(2) : '1.00'; }
				update();
			});

			// Fetch TSV and initialize
			fetch(TSV_PATH).then(r => {
				if (!r.ok) throw new Error('Failed to fetch TSV: ' + r.status);
				return r.text();
			}).then(text => {
				const parsed = parseTSV(text);
				header = parsed.header;
				rows = parsed.data;
				// initial population: set full option lists
				LEVELS.forEach(l => {
					const set = new Set(rows.map(r => r[l]).filter(Boolean));
					const arr = Array.from(set).sort();
					// populate datalist for the search input
					setSelectOptions(selects[l], arr);
				});
				// ensure child options reflect any top-level defaults (none at start)
				refreshAllOptions();
				update();
			}).catch(err => {
				document.getElementById('plot').textContent = 'Error loading data: ' + err.message;
				console.error(err);
			});

			// Build TSV string from filtered rows using header order
			function buildTSVFromRows(filteredRows) {
				if (!header || header.length === 0) {
					// fall back to keys from first row
					if (filteredRows.length === 0) return '';
					header = Object.keys(filteredRows[0]);
				}
				const lines = [];
				lines.push(header.join('\t'));
				filteredRows.forEach(r => {
					const vals = header.map(h => {
						let v = r[h];
						if (v === undefined || v === null) return '';
						// escape tabs/newlines
						return String(v).replace(/\t/g, ' ').replace(/\n/g, ' ');
					});
					lines.push(vals.join('\t'));
				});
				return lines.join('\n');
			}

			// Trigger a download of current filtered TSV
			function downloadFilteredTSV() {
				const filters = getFilters();
				const filtered = filterData(rows, filters);
				const tsv = buildTSVFromRows(filtered);
				if (!tsv) {
					alert('No data to download');
					return;
				}
				const blob = new Blob([tsv], { type: 'text/tab-separated-values;charset=utf-8;' });
				const url = URL.createObjectURL(blob);
				const a = document.createElement('a');
				const now = new Date().toISOString().slice(0,19).replace(/[:T]/g,'-');
				a.href = url;
				a.download = `arsc_filtered_${now}.tsv`;
				document.body.appendChild(a);
				a.click();
				document.body.removeChild(a);
				URL.revokeObjectURL(url);
			}

			if (downloadBtn) downloadBtn.addEventListener('click', downloadFilteredTSV);
