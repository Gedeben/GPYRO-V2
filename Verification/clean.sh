#!/bin/bash

BASE_DIR="./"
LONG_CASES=("0.0001mm_dt0.05s" "0.001mm_dt0.05s" "0.1mm_dt0.0005s")
OMP_FOLDER="OMP_parallelization"


echo "🔎 Scanning and cleaning in: $BASE_DIR"

# Extensions et fichiers autorisés
allowed_exts=("sh" "py" "data" "tex" "bib" "txt" "pdf" "fds" "md" "cmp" "cnd" "ini" "ssf")
allowed_csv=("analytical_results.csv" "heat_conduction_kc.csv" "reference_cc_devc.csv" "output_ThermaKin.csv" "1d_rad_conv_devc.csv" "radiation_loss_devc.csv" "MaCFP-PMMA_Gasification_q50_Mass_R3.csv")

# Liste des fichiers à supprimer
to_delete=()

while IFS= read -r -d '' file; do
  filename=$(basename "$file")
  extension="${filename##*.}"
  dirpath=$(dirname "$file")
  keep=false

  # 🔒 Ne jamais supprimer les *.summary.csv ou *timing.csv dans les cas longs
  for long_dir in "${LONG_CASES[@]}"; do
    if [[ "$dirpath" == *"$long_dir"* ]]; then
      if [[ "$filename" == *summary*.csv || "$filename" == *timing.csv ]]; then
        keep=true
        break
      fi
    fi
  done
  if $keep; then
    continue
  fi

  # ✅ Garde les extensions autorisées
  for ext in "${allowed_exts[@]}"; do
    if [[ "$filename" == *.$ext ]]; then
      keep=true
      break
    fi
  done

  # ✅ Garde les CSV explicitement autorisés
  if [[ "$extension" == "csv" && $keep == false ]]; then
    for allowed_name in "${allowed_csv[@]}"; do
      if [[ "$filename" == "$allowed_name" ]]; then
        keep=true
        break
      fi
    done
  fi

  if ! $keep; then
    to_delete+=("$file")
  fi
done < <(find "$BASE_DIR" -type f -print0)


# 🔍 Chercher le dossier OMP_FOLDER (partiel match) et lister ses sous-dossiers
omp_target_dirs=()
while IFS= read -r -d '' match_dir; do
  while IFS= read -r -d '' subdir; do
    omp_target_dirs+=("$subdir")
  done < <(find "$match_dir" -mindepth 1 -maxdepth 1 -type d -print0)
done < <(find "$BASE_DIR" -type d -name "*$OMP_FOLDER*" -print0)

if [[ ${#omp_target_dirs[@]} -gt 0 ]]; then
  echo "🗑️ Folders inside '$OMP_FOLDER' to be deleted:"
  for d in "${omp_target_dirs[@]}"; do
    echo "  $d"
  done
fi



# 🗑️ Liste des fichiers à supprimer
echo "🗑️ Files to be deleted:"
for f in "${to_delete[@]}"; do
  echo "  $f"
done

# Demande confirmation
read -p "❓ Do you want to delete these files? [y/N] " confirm
if [[ "$confirm" == [yY] ]]; then
  for f in "${to_delete[@]}"; do
    rm -f "$f"
  done
  for d in "${omp_target_dirs[@]}"; do
    rm -rf "$d"
  done
  echo "✅ Files deleted."
else
  echo "❌ Deletion canceled."
fi



