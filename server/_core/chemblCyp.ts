import { z } from "zod";

// ChEMBL public REST API base
const CHEMBL_BASE = "https://www.ebi.ac.uk/chembl/api/data";

export type CypIsoform =
  | "CYP1A2"
  | "CYP2A6"
  | "CYP2B6"
  | "CYP2C8"
  | "CYP2C9"
  | "CYP2C19"
  | "CYP2D6"
  | "CYP2E1"
  | "CYP3A4"
  | "CYP3A5";

export const CYP_TARGETS: Record<CypIsoform, string> = {
  CYP1A2: "CHEMBL3356",
  CYP2A6: "CHEMBL5282",
  CYP2B6: "CHEMBL4729",
  CYP2C8: "CHEMBL3721",
  CYP2C9: "CHEMBL3397",
  CYP2C19: "CHEMBL3622",
  CYP2D6: "CHEMBL289",
  CYP2E1: "CHEMBL5281",
  CYP3A4: "CHEMBL340",
  CYP3A5: "CHEMBL3019",
};

const ChemblActivitySchema = z.object({
  activity_id: z.number().optional(),
  molecule_chembl_id: z.string().nullable().optional(),
  target_chembl_id: z.string().nullable().optional(),
  assay_chembl_id: z.string().nullable().optional(),
  document_chembl_id: z.string().nullable().optional(),
  standard_type: z.string().nullable().optional(),
  standard_relation: z.string().nullable().optional(),
  standard_value: z.string().nullable().optional(),
  standard_units: z.string().nullable().optional(),
  pchembl_value: z.number().nullable().optional(),
});

const ChemblActivityRespSchema = z.object({
  page_meta: z.any().optional(),
  activities: z.array(ChemblActivitySchema).optional(),
});

export interface CypEvidenceBest {
  isoform: string;
  targetChEMBLId: string;
  best: {
    standardType: string;
    standardValue: number;
    standardUnits: string;
    relation: string | null;
    pchemblValue: number | null;
    assayChEMBLId: string | null;
    documentChEMBLId: string | null;
    moleculeChEMBLId: string | null;
  } | null;
  evidenceCount: number;
  samples: Array<{
    standardType: string;
    standardValue: number;
    standardUnits: string;
    relation: string | null;
    pchemblValue: number | null;
    assayChEMBLId: string | null;
    documentChEMBLId: string | null;
    moleculeChEMBLId: string | null;
  }>;
}

function normalizeValueToNm(value: number, units: string): number | null {
  const u = (units || "").toLowerCase();
  if (u === "nm") return value;
  if (u === "um" || u === "µm") return value * 1000;
  if (u === "pm") return value / 1000;
  // Sometimes ChEMBL uses 'M' or 'nM' like strings; we handle above.
  return null;
}

function pickBestActivity(rows: Array<{ valueNm: number; raw: any }>) {
  if (rows.length === 0) return null;
  // Best = smallest potency value (IC50/Ki) in nM
  rows.sort((a, b) => a.valueNm - b.valueNm);
  return rows[0];
}

async function fetchJson(url: string) {
  const res = await fetch(url, {
    headers: { Accept: "application/json" },
    signal: AbortSignal.timeout(20000),
  });
  if (!res.ok) {
    throw new Error(`ChEMBL HTTP ${res.status}`);
  }
  return res.json();
}

export async function fetchChemblCypEvidenceByMoleculeChemblId(opts: {
  moleculeChemblId: string;
  isoforms?: CypIsoform[];
  /** Total activities pulled per type (IC50 / Ki). */
  limitTotal?: number;
}): Promise<CypEvidenceBest[]> {
  const isoforms = opts.isoforms ?? (Object.keys(CYP_TARGETS) as CypIsoform[]);
  const limit = Math.min(Math.max(opts.limitTotal ?? 200, 50), 500);

  // Evidence-first strategy (efficient):
  // Pull activities for the molecule (IC50 + Ki), then filter by CYP target ids.
  const urls = [
    `${CHEMBL_BASE}/activity?molecule_chembl_id=${encodeURIComponent(
      opts.moleculeChemblId
    )}&standard_type=IC50&limit=${limit}`,
    `${CHEMBL_BASE}/activity?molecule_chembl_id=${encodeURIComponent(
      opts.moleculeChemblId
    )}&standard_type=Ki&limit=${limit}`,
  ];

  const allRows: any[] = [];
  for (const u of urls) {
    try {
      const data = ChemblActivityRespSchema.parse(await fetchJson(u));
      if (Array.isArray((data as any).activities)) {
        allRows.push(...((data as any).activities as any[]));
      }
    } catch {
      // ignore
    }
  }

  const results: CypEvidenceBest[] = [];

  for (const iso of isoforms) {
    const target = CYP_TARGETS[iso];
    const rows = allRows.filter(
      r => String(r.target_chembl_id ?? "") === target
    );

    const parsed = rows
      .map(r => {
        const standardType = String(r.standard_type ?? "").trim();
        const standardUnits = String(r.standard_units ?? "").trim();
        const rel = r.standard_relation ? String(r.standard_relation) : null;
        const v = r.standard_value != null ? Number(r.standard_value) : NaN;
        if (!Number.isFinite(v)) return null;
        const nm = normalizeValueToNm(v, standardUnits);
        if (nm == null) return null;
        return {
          valueNm: nm,
          raw: {
            standardType,
            standardUnits,
            standardValue: v,
            relation: rel,
            pchemblValue:
              r.pchembl_value != null &&
              Number.isFinite(Number(r.pchembl_value))
                ? Number(r.pchembl_value)
                : null,
            assayChEMBLId: r.assay_chembl_id ?? null,
            documentChEMBLId: r.document_chembl_id ?? null,
            moleculeChEMBLId: r.molecule_chembl_id ?? null,
          },
        };
      })
      .filter(Boolean) as Array<{ valueNm: number; raw: any }>;

    const best = pickBestActivity(parsed);
    const samples = parsed
      .sort((a, b) => a.valueNm - b.valueNm)
      .slice(0, 5)
      .map(x => x.raw);

    results.push({
      isoform: iso,
      targetChEMBLId: target,
      best: best ? { ...best.raw } : null,
      evidenceCount: parsed.length,
      samples,
    });
  }

  return results;
}

export async function resolveMoleculeChemblIdBySmiles(smiles: string) {
  const s = smiles.trim();
  if (!s) return null;

  // Prefer similarity query (more robust than exact), keep first hit.
  // docs: /similarity/{smiles}/{similarity} or /similarity/{smiles}
  const url = `${CHEMBL_BASE}/similarity/${encodeURIComponent(
    s
  )}?similarity=100&limit=1`;

  const data = await fetchJson(url);
  const molecules = data?.molecules;
  const id = molecules?.[0]?.molecule_chembl_id ?? null;
  return typeof id === "string" ? id : null;
}
