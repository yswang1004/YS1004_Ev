import { useMemo, useState } from "react";
import { Navbar } from "@/components/Navbar";
import { Card, CardContent } from "@/components/ui/card";
import { Button } from "@/components/ui/button";
import { Input } from "@/components/ui/input";
import { Alert } from "@/components/ui/alert";

type CalibRow = {
  compound_name: string;
  ic50_uM: number;
  reference?: string;
  assay?: string;
};

function parseCSV(text: string): CalibRow[] {
  const lines = text
    .split(/\r?\n/)
    .map(l => l.trim())
    .filter(Boolean);
  if (lines.length < 2) return [];

  const header = lines[0].split(",").map(h => h.trim().toLowerCase());
  const idx = (name: string) => header.indexOf(name);

  const iName = idx("compound_name");
  const iIc50 = idx("ic50_um");
  if (iName < 0 || iIc50 < 0) return [];

  const opt = {
    reference: idx("reference"),
    assay: idx("assay"),
  };

  const rows: CalibRow[] = [];
  for (let i = 1; i < lines.length; i++) {
    const cols = lines[i].split(",");
    const name = (cols[iName] ?? "").trim();
    if (!name) continue;
    const raw = (cols[iIc50] ?? "").trim();
    const ic50 = Number(raw);
    if (!Number.isFinite(ic50)) continue;

    const row: CalibRow = { compound_name: name, ic50_uM: ic50 };

    (Object.keys(opt) as (keyof typeof opt)[]).forEach(k => {
      const j = opt[k];
      if (j != null && j >= 0) {
        const v = (cols[j] ?? "").trim();
        if (v) (row as any)[k] = v;
      }
    });

    rows.push(row);
  }
  return rows;
}

export default function CYP2E1Calibration() {
  const [raw, setRaw] = useState<string>(
    "compound_name,ic50_uM,assay,reference\nExampleCompound,1.2,,"
  );
  const [error, setError] = useState<string>("");
  const [saved, setSaved] = useState<boolean>(false);

  const parsed = useMemo(() => {
    try {
      const rows = parseCSV(raw);
      setError("");
      return rows;
    } catch {
      setError("Failed to parse CSV. Please check formatting.");
      return [];
    }
  }, [raw]);

  const onFile = async (file: File) => {
    const text = await file.text();
    setRaw(text);
    setSaved(false);
  };

  const save = () => {
    if (parsed.length === 0) {
      setError("No valid rows found. CSV must include: compound_name, ic50_uM");
      return;
    }
    localStorage.setItem("cyp2e1Calibration", JSON.stringify(parsed));
    setSaved(true);
    setError("");
  };

  const clear = () => {
    localStorage.removeItem("cyp2e1Calibration");
    setSaved(false);
  };

  return (
    <div className="min-h-screen bg-background text-foreground">
      <Navbar />
      <main className="container py-8 space-y-6">
        <div className="space-y-1">
          <h1 className="text-2xl font-semibold tracking-tight">
            CYP2E1 Calibration (IC50, µM)
          </h1>
          <p className="text-sm text-muted-foreground">
            Upload experimental IC50 data (µM) to reduce false positives. Results
            will adjust rank automatically.
          </p>
        </div>

        {error && (
          <Alert className="border-rose-500/40 bg-rose-500/10 text-rose-200">
            {error}
          </Alert>
        )}
        {saved && (
          <Alert className="border-emerald-500/40 bg-emerald-500/10 text-emerald-200">
            Saved. Results page will apply calibration automatically.
          </Alert>
        )}

        <Card className="card-glow bg-card border-border">
          <CardContent className="p-6 space-y-4">
            <div className="flex items-center gap-3">
              <Input
                type="file"
                accept=".csv,text/csv"
                onChange={e => {
                  const f = e.target.files?.[0];
                  if (f) onFile(f);
                }}
              />
              <Button onClick={save}>Save</Button>
              <Button variant="secondary" onClick={clear}>
                Clear
              </Button>
            </div>

            <div className="space-y-2">
              <div className="text-sm font-medium">CSV Preview</div>
              <textarea
                className="w-full h-64 rounded-md border border-border bg-background p-3 font-mono text-xs"
                value={raw}
                onChange={e => {
                  setRaw(e.target.value);
                  setSaved(false);
                }}
              />
              <div className="text-xs text-muted-foreground">
                Required columns: <code>compound_name</code>, <code>ic50_uM</code>.
                Optional: <code>assay,reference</code>.
              </div>
            </div>

            <div className="text-sm text-muted-foreground">
              Parsed rows: <span className="text-foreground">{parsed.length}</span>
            </div>
          </CardContent>
        </Card>
      </main>
    </div>
  );
}
