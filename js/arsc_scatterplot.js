// Path to the TSV file (relative)
const TSV_PATH = './data/arsc_gtdb_r226rep.tsv';

// Hierarchy levels (top -> bottom)
const LEVELS = ['domain','phylum','class','order','family','genus'];

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

// loading indicator element (shown during initial TSV fetch)
const loadingEl = document.getElementById('loading');
// file input for user FASTA
let fastaFileInput = document.getElementById('fastaFile');
let fastaNameSpan = document.getElementById('fastaName');
const fastaContainer = document.getElementById('fastaFile') ? document.getElementById('fastaFile').closest('div') : null;
const origFastaContainerHTML = fastaContainer ? fastaContainer.innerHTML : null;
let parseWorker = null;
let userSample = null; // stores last computed FASTA result { filename, N_ARSC, C_ARSC, S_ARSC }

function bindFastaElements() {
	// re-query elements (useful after replacing container innerHTML)
	fastaFileInput = document.getElementById('fastaFile');
	fastaNameSpan = document.getElementById('fastaName');
	if (!fastaFileInput) return;
	// remove previous listeners by cloning
	const old = fastaFileInput;
	const newInput = old.cloneNode(true);
	old.parentNode.replaceChild(newInput, old);
	fastaFileInput = newInput;
	fastaFileInput.addEventListener('change', handleFastaFileChange);
}

function showUserSampleInfo() {
	if (!fastaContainer || !userSample) return;
	let txt = `file: ${userSample.filename}; N-ARSC: ${userSample.N_ARSC.toFixed(5)}; C-ARSC: ${userSample.C_ARSC.toFixed(5)}; S-ARSC: ${userSample.S_ARSC.toFixed(5)}`;
	let html = `<div style="font-size:0.95rem;color:#222;">${txt}</div>`;
	if (userSample.warning) {
		html += `<div style="color:#b33;margin-top:6px;font-size:0.85rem;">⚠ ${userSample.warning.message}</div>`;
	}
	fastaContainer.innerHTML = html;
}

function restoreFastaContainer() {
	if (!fastaContainer || !origFastaContainerHTML) return;
	fastaContainer.innerHTML = origFastaContainerHTML;
	bindFastaElements();
}

// helper to add/remove previous user-sample overlay
function removeUserSampleOverlay() {
	// remove traces named 'user-sample' and annotations with user-sample id
	try {
		const gd = document.getElementById('plot');
		if (!gd || !gd.data) return;
		const tracesToRemove = [];
		gd.data.forEach((t, idx) => { if (t && t.name === 'user-sample') tracesToRemove.push(idx); });
		// remove from last to first to keep indices valid
		for (let i = tracesToRemove.length - 1; i >= 0; i--) {
			Plotly.deleteTraces(gd, tracesToRemove[i]);
		}
		// remove annotations added with uid 'user-sample-annotation'
		const layout = gd.layout || {};
		if (layout.annotations) {
			const kept = layout.annotations.filter(a => !(a && a._userSample));
			Plotly.relayout(gd, { annotations: kept });
		}
			// remove shapes (lines) added for user-sample
			if (layout.shapes) {
				const keptShapes = layout.shapes.filter(s => !(s && s._userSample));
				Plotly.relayout(gd, { shapes: keptShapes });
			}
	} catch (e) { console.warn(e); }
}

// add overlay (horizontal dashed line + annotation) for currently stored userSample
function addUserSampleOverlay() {
	if (!userSample) return;
	try {
		const gd = document.getElementById('plot');
		if (!gd) return;
		const yFieldCur = (ySelect && ySelect.value) ? ySelect.value : 'N_ARSC';
		const yValCur = (userSample && userSample[yFieldCur] !== undefined) ? userSample[yFieldCur] : userSample.N_ARSC;
		const existingShapes = (gd.layout && gd.layout.shapes) ? gd.layout.shapes.slice() : [];
		const line = { type: 'line', xref: 'paper', x0: 0, x1: 1, yref: 'y', y0: yValCur, y1: yValCur, line: { color: 'red', width: 2 }, _userSample: true };
		Plotly.relayout('plot', { shapes: existingShapes.concat([line]) });
		const existingAnns = (gd.layout && gd.layout.annotations) ? gd.layout.annotations.slice() : [];
		const ann = { xref: 'paper', x: 0.99, xanchor: 'right', yref: 'y', y: yValCur, text: `${userSample.filename} — ${yFieldCur}: ${yValCur.toFixed(5)}`, bgcolor: '#fff8', bordercolor: '#ff3333', font: { color: '#800', size: 12 }, _userSample: true };
		Plotly.relayout('plot', { annotations: existingAnns.concat([ann]) });
	} catch (e) { console.warn(e); }
}

// handle file selection
function handleFastaFileChange(ev) {
	const f = ev.target.files && ev.target.files[0];
	if (!f) return;

	// validate extension early to avoid unnecessary parsing
	const allowedExts = ['.fa', '.faa', '.fasta'];
	const nameLower = (f.name || '').toLowerCase();
	if (!allowedExts.some(ext => nameLower.endsWith(ext))) {
		try { alert('Unsupported file type. Please upload a .fa, .faa or .fasta file.'); } catch (e) {}
		// clear input so user can try again
		if (fastaFileInput) { try { fastaFileInput.value = ''; } catch (e) {} }
		if (fastaNameSpan) fastaNameSpan.textContent = '';
		return;
	}
	if (fastaNameSpan) fastaNameSpan.textContent = `${f.name} (${(f.size/1024).toFixed(1)} KB)`;
	const MAX_BYTES = 10 * 1024 * 1024; // 10MB
	if (f.size > MAX_BYTES) {
		const ok = confirm('File seems large (>10MB). Continue parsing in browser? Recommended: use smaller file. Continue?');
		if (!ok) return;
	}
	// read as text and spawn worker
	const reader = new FileReader();
	reader.onload = function(e) {
		const text = e.target.result;
		// start worker
		if (parseWorker) {
			parseWorker.terminate();
			parseWorker = null;
		}
		parseWorker = new Worker('js/arsc_worker.js');
		if (loadingEl) { loadingEl.style.display = 'inline-block'; loadingEl.textContent = 'Parsing FASTA...'; }
		parseWorker.onmessage = function(me) {
			const m = me.data;
			if (!m || !m.type) return;
			if (m.type === 'result') {
				if (loadingEl) loadingEl.style.display = 'none';
				// Handle warnings from worker: strong check for nucleotide-like input
				if (m.warning && m.warning.type === 'nucleotide_detected') {
					// prompt user: proceed or cancel
					let proceed = false;
					try {
						proceed = confirm(
							"⚠️ Warning: Uploaded sequence looks like DNA/RNA (mostly A/T/G/C/U/N). Results may be incorrect.\n" +
							"Do you want to continue? (OK = continue, Cancel = cancel)"
						);
					} catch (e) { proceed = false; }
					if (!proceed) {
						// user cancelled: stop processing and reset file input
						if (loadingEl) loadingEl.style.display = 'none';
						try { parseWorker.terminate(); } catch (e) {}
						parseWorker = null;
						if (fastaFileInput) { try { fastaFileInput.value = ''; } catch (e) {} }
						if (fastaNameSpan) fastaNameSpan.textContent = '';
						// leave userSample unchanged (do not store)
						return;
					}
				}
				// store the computed sample so it persists across plot redraws
				userSample = {
					filename: m.filename || 'uploaded',
					N_ARSC: typeof m.N_ARSC === 'number' ? m.N_ARSC : parseFloat(m.N_ARSC),
					C_ARSC: typeof m.C_ARSC === 'number' ? m.C_ARSC : parseFloat(m.C_ARSC),
					S_ARSC: typeof m.S_ARSC === 'number' ? m.S_ARSC : parseFloat(m.S_ARSC),
					warning: m.warning || null
				};
				// for non-nucleotide warnings (unknown residues), show a note
				if (m.warning && m.warning.type === 'unknown_residues') {
					try { alert(`Warning: ${m.warning.message}`); } catch (e) {}
				}
				// replace the input UI with summary info
				showUserSampleInfo();
				// update the plot (drawPlot will re-apply overlays)
				try { update(); } catch (err) { /* ignore */ }
				// terminate worker
				parseWorker.terminate(); parseWorker = null;
			} else if (m.type === 'error') {
				if (loadingEl) loadingEl.style.display = 'none';
				alert('Error parsing FASTA: ' + m.message);
				parseWorker.terminate(); parseWorker = null;
			}
		};
		parseWorker.postMessage({ type: 'parse', text, filename: f.name });
	};
	reader.onerror = function(err) { alert('Failed to read file: ' + err); };
	reader.readAsText(f);
}

// bind elements initially
bindFastaElements();

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
		
		// Re-attach click handler after plot recreation
		setupPlotClickHandler();
		
		// if a user sample exists, re-apply its overlay after the main plot redraw
		if (userSample) addUserSampleOverlay();
	}// Show modal with copyable data
function showDataModal(row, xField, yField) {
	// Build formatted text content
	const taxLevels = ['domain', 'phylum', 'class', 'order', 'family', 'genus'];
	let text = '';
	
	if (row['id']) text += `ID: ${row['id']}\n`;
	text += `\nTaxonomy:\n`;
	taxLevels.forEach(level => {
		if (row[level]) text += `  ${level}: ${row[level]}\n`;
	});
	
	text += `\nARSC Values:\n`;
	text += `  N-ARSC: ${row['N_ARSC']}\n`;
	text += `  C-ARSC: ${row['C_ARSC']}\n`;
	text += `  S-ARSC: ${row['S_ARSC']}\n`;
	
	text += `\nOther:\n`;
	text += `  ${xField}: ${row[xField]}\n`;
	if (row['AvgResMW']) text += `  AvgResMW: ${row['AvgResMW']}\n`;
	if (row['sum_len']) text += `  sum_len: ${row['sum_len']}\n`;
	
	// Create modal
	const modal = document.createElement('div');
	modal.style.cssText = `
		position: fixed;
		top: 0;
		left: 0;
		width: 100%;
		height: 100%;
		background: rgba(0,0,0,0.5);
		display: flex;
		align-items: center;
		justify-content: center;
		z-index: 10000;
	`;
	
	const modalContent = document.createElement('div');
	modalContent.style.cssText = `
		background: white;
		padding: 20px;
		border-radius: 8px;
		max-width: 500px;
		width: 90%;
		max-height: 80vh;
		overflow: auto;
		box-shadow: 0 4px 6px rgba(0,0,0,0.1);
	`;
	
	const title = document.createElement('h3');
	title.textContent = 'Sample Data';
	title.style.marginTop = '0';
	
	const textarea = document.createElement('textarea');
	textarea.value = text;
	textarea.style.cssText = `
		width: 100%;
		height: 300px;
		font-family: monospace;
		font-size: 13px;
		padding: 10px;
		border: 1px solid #ccc;
		border-radius: 4px;
		resize: vertical;
		box-sizing: border-box;
	`;
	textarea.readOnly = true;
	
	const buttonContainer = document.createElement('div');
	buttonContainer.style.cssText = `
		margin-top: 15px;
		display: flex;
		gap: 10px;
		justify-content: flex-end;
	`;
	
	const copyBtn = document.createElement('button');
	copyBtn.textContent = 'Copy';
	copyBtn.style.cssText = `
		padding: 8px 16px;
		background: #007bff;
		color: white;
		border: none;
		border-radius: 4px;
		cursor: pointer;
	`;
	copyBtn.onclick = () => {
		textarea.select();
		document.execCommand('copy');
		copyBtn.textContent = 'Copied!';
		setTimeout(() => { copyBtn.textContent = 'Copy'; }, 1500);
	};
	
	const closeBtn = document.createElement('button');
	closeBtn.textContent = 'Close';
	closeBtn.style.cssText = `
		padding: 8px 16px;
		background: #6c757d;
		color: white;
		border: none;
		border-radius: 4px;
		cursor: pointer;
	`;
	closeBtn.onclick = () => {
		document.body.removeChild(modal);
	};
	
	buttonContainer.appendChild(copyBtn);
	buttonContainer.appendChild(closeBtn);
	
	modalContent.appendChild(title);
	modalContent.appendChild(textarea);
	modalContent.appendChild(buttonContainer);
	modal.appendChild(modalContent);
	
	// Close on background click
	modal.onclick = (e) => {
		if (e.target === modal) {
			document.body.removeChild(modal);
		}
	};
	
	// Close on Escape key
	const escHandler = (e) => {
		if (e.key === 'Escape') {
			document.body.removeChild(modal);
			document.removeEventListener('keydown', escHandler);
		}
	};
	document.addEventListener('keydown', escHandler);
	
	document.body.appendChild(modal);
	
	// Auto-select text for easy copying
	setTimeout(() => { textarea.select(); }, 100);
}

// Update routine: get filters, filter rows, draw
function update() {
	const filters = getFilters();
	const filtered = filterData(rows, filters);
	drawPlot(filtered);
}

// Setup click event listener once (to avoid memory leak from repeated registration)
function setupPlotClickHandler() {
	const plotDiv = document.getElementById('plot');
	if (!plotDiv) return;
	
	plotDiv.on('plotly_click', function(data) {
		if (!data.points || data.points.length === 0) return;
		const pt = data.points[0];
		
		// Get current axis fields
		const yField = ySelect && ySelect.value ? ySelect.value : 'N_ARSC';
		const xField = 'GC(%)';
		
		// Get the original data row for this point
		const filters = getFilters();
		const filteredRows = filterData(rows, filters);
		
		// Find the corresponding row by matching x, y values
		const xVal = pt.x;
		const yVal = pt.y;
		const matchedRow = filteredRows.find(r => 
			Math.abs(r[xField] - xVal) < 0.0001 && 
			Math.abs(r[yField] - yVal) < 0.0001
		);
		
		if (matchedRow) {
			showDataModal(matchedRow, xField, yField);
		}
	});
}

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
	const listId = sel.id.endsWith('Search') ? sel.id.replace('Search', 'List') : sel.id.replace('Select', 'List');
	const dl = document.getElementById(listId);
	if (dl) {
		while (dl.firstChild) dl.removeChild(dl.firstChild);
		values.forEach(v => { const o = document.createElement('option'); o.value = v; dl.appendChild(o); });
	}
	// restore if still available
	if (cur && values.includes(cur)) sel.value = cur;
	else sel.value = '';
}

// Refresh all select options based on current ancestor selections
function refreshAllOptions() {
	LEVELS.forEach((level) => {
		const vals = availableOptionsForLevel(level);
		setSelectOptions(selects[level], vals);
		const s = selects[level];
		if (s && s.value && !vals.includes(s.value)) s.value = '';
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
	setParentsFromChild(level);
	// After setting parents, refresh options to reflect consistent children sets
	refreshAllOptions();
	update();
}

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

// ========== Initialization ==========

// assign data-level attributes for easier handling
LEVELS.forEach(l => { selects[l].dataset.level = l; });

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
		} else {
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
	if (markerSizeInput) markerSizeInput.value = 8;
	if (markerSizeVal) markerSizeVal.textContent = '8';
	if (markerAlphaInput) markerAlphaInput.value = 1.0;
	if (markerAlphaVal) markerAlphaVal.textContent = '1.00';
	// clear stored user sample and input UI as part of reset
	userSample = null;
	removeUserSampleOverlay();
	// restore the original upload UI and clear any displayed filename/text
	restoreFastaContainer();
	if (fastaFileInput) fastaFileInput.value = '';
	if (fastaNameSpan) fastaNameSpan.textContent = '';
	update();
});

if (downloadBtn) downloadBtn.addEventListener('click', downloadFilteredTSV);

// Fetch TSV and initialize
try {
	if (loadingEl) loadingEl.style.display = 'inline-block';
} catch (e) { /* ignore */ }

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
	// Setup click handler once after initial plot
	setupPlotClickHandler();
	if (loadingEl) loadingEl.style.display = 'none';
}).catch(err => {
	if (loadingEl) loadingEl.style.display = 'none';
	document.getElementById('plot').textContent = 'Error loading data: ' + err.message;
	console.error(err);
});
