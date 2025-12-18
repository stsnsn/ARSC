// Worker to parse FASTA text and compute N/C/S-ARSC per residue (fast version)
// Receives: { type: 'parse', text: string, filename: string }
// Posts messages: { type: 'result', filename, N_ARSC, C_ARSC, S_ARSC, warning? } or { type: 'error', message }

// --- 高速化のためのルックアップテーブル（ASCII 0-127） ---
const TABLE_N = new Int8Array(128).fill(-1);
const TABLE_C = new Int8Array(128).fill(-1);
const TABLE_S = new Int8Array(128).fill(-1);
const TABLE_IS_NUC = new Int8Array(128).fill(0); // table to quickly test whether a character is a nucleotide (A,T,G,C,U,N)

const RESIDUE_DATA = {
  A: { N:0, C:1, S:0 }, R: { N:3, C:4, S:0 }, N: { N:1, C:2, S:0 },
  D: { N:0, C:2, S:0 }, C: { N:0, C:1, S:1 }, Q: { N:1, C:3, S:0 },
  E: { N:0, C:3, S:0 }, G: { N:0, C:0, S:0 }, H: { N:2, C:4, S:0 },
  I: { N:0, C:4, S:0 }, L: { N:0, C:4, S:0 }, K: { N:1, C:4, S:0 },
  M: { N:0, C:3, S:1 }, F: { N:0, C:7, S:0 }, P: { N:0, C:3, S:0 },
  S: { N:0, C:1, S:0 }, T: { N:0, C:2, S:0 }, W: { N:1, C:9, S:0 },
  Y: { N:0, C:7, S:0 }, V: { N:0, C:3, S:0 }
};

// initialize tables for both upper and lower case ASCII
for (const [aa, counts] of Object.entries(RESIDUE_DATA)) {
  const up = aa.charCodeAt(0);
  const lo = aa.toLowerCase().charCodeAt(0);
  TABLE_N[up] = counts.N; TABLE_N[lo] = counts.N;
  TABLE_C[up] = counts.C; TABLE_C[lo] = counts.C;
  TABLE_S[up] = counts.S; TABLE_S[lo] = counts.S;
}

// define nucleotide characters (A,T,G,C,U,N) for upper and lower case
['A','T','G','C','U','N'].forEach(ch => {
  const up = ch.charCodeAt(0);
  const lo = ch.toLowerCase().charCodeAt(0);
  TABLE_IS_NUC[up] = 1;
  TABLE_IS_NUC[lo] = 1;
});

function parseFASTAandComputeFast(text) {
  let totalResidues = 0;
  let N_total = 0, C_total = 0, S_total = 0;
  let unknownCount = 0;
  // nucleotide detection counters
  let validCharCount = 0; // count of header-excluded valid characters considered
  let nucleotideMatchCount = 0; // count of characters matching A/T/G/C/U/N

  // 不明な文字の種類を記録（ASCII範囲）
  const unknownFlags = new Uint8Array(128);

  const len = text.length;
  let isHeader = false;

  for (let i = 0; i < len; i++) {
    const code = text.charCodeAt(i);

    // 改行処理
    if (code === 10 || code === 13) { isHeader = false; continue; }
    if (isHeader) continue;
    if (code === 62) { isHeader = true; continue; } // '>'

    // ASCII 範囲外の文字コードは無視
    if (code > 127) continue;
    // 制御文字や空白、数字なども無視
    if (code < 65) continue;

    // count this as a valid (non-control, non-header) character for nucleotide detection
    validCharCount++;
    if (TABLE_IS_NUC[code] === 1) nucleotideMatchCount++;

    const n = TABLE_N[code];
    if (n !== -1) {
      totalResidues++;
      N_total += n;
      C_total += TABLE_C[code];
      S_total += TABLE_S[code];
    } else {
      // 定義されていない文字（かつ code >= 65）の場合、アルファベットのみを警告対象にする
      if ((code >= 65 && code <= 90) || (code >= 97 && code <= 122)) {
        unknownCount++;
        unknownFlags[code] = 1;
      }
    }
  }

  if (totalResidues === 0) {
    throw new Error('No valid residues found.');
  }

  let warning = null;
  // nucleotide-detection: if >=98% of valid characters are A/T/G/C/U/N, treat as likely nucleotide sequence
  if (validCharCount > 0) {
    const ratio = nucleotideMatchCount / validCharCount;
    if (ratio > 0.98) {
      warning = {
        type: 'nucleotide_detected',
        message: 'Caution: The input appears to be DNA/RNA (mostly A,T,G,C,U or N). This tool expects amino-acid sequences; results may be meaningless.'
      };
    }
  }
  // if not flagged as nucleotide, fall back to unknown-residue warning
  if (!warning && unknownCount > 0) {
    const foundChars = [];
    for (let c = 0; c < 128; c++) if (unknownFlags[c] === 1) foundChars.push(String.fromCharCode(c));
    warning = { type: 'unknown_residues', count: unknownCount, chars: foundChars, message: `Excluded ${unknownCount} unknown residue(s): ${foundChars.join(', ')}` };
  }

  return {
    N_ARSC: N_total / totalResidues,
    C_ARSC: C_total / totalResidues,
    S_ARSC: S_total / totalResidues,
    totalResidues,
    warning
  };
}

onmessage = function(ev) {
  const msg = ev.data;
  if (!msg || msg.type !== 'parse') return;
  try {
    const text = msg.text || '';
    const res = parseFASTAandComputeFast(text);
    postMessage({ type: 'result', filename: msg.filename || 'uploaded', N_ARSC: res.N_ARSC, C_ARSC: res.C_ARSC, S_ARSC: res.S_ARSC, warning: res.warning });
  } catch (err) {
    postMessage({ type: 'error', message: err.message || String(err) });
  }
};
