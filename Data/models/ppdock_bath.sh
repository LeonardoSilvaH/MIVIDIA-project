#!/bin/bash

# ===== CONFIGURAÇÕES (ajuste se necessário) =====
PPDOCK_REPO="/hd/leo/PPDOCK"                        # clone do GitHub
CONDA_ENV="PPDOCK"
TEMP_BASE="/tmp/ppdock_work"                       # base para trabalhos temporários
RESULT_DIR="/hd/leo/models/resultados_docking"     # onde salvar as saídas finais
# =================================================

# Ativa o ambiente conda
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV"

# Cria a pasta de resultados se não existir
mkdir -p "$RESULT_DIR"

# Percorre todas as pastas de proteína (primeiro nível)
for protein_dir in */; do
    protein_dir="${protein_dir%/}"
    [ ! -d "$protein_dir" ] && continue

    # Percorre as mutações (subpastas dentro da proteína)
    for mut_dir in "$protein_dir"/*/; do
        mut_dir="${mut_dir%/}"
        [ ! -d "$mut_dir" ] && continue

        # Localiza o primeiro arquivo .pdb dentro da mutação
        pdb_file=$(ls "$mut_dir"/*.pdb 2>/dev/null | head -n 1)
        if [ -z "$pdb_file" ]; then
            echo "Aviso: Nenhum arquivo .pdb encontrado em $mut_dir, pulando..."
            continue
        fi

        mutation_name=$(basename "$mut_dir")

        # Lista todos os ligantes .sdf na raiz da proteína
        for ligand_sdf in "$protein_dir"/*.sdf; do
            [ -f "$ligand_sdf" ] || continue
            ligand_name=$(basename "$ligand_sdf" .sdf)

            echo "Processando: Proteina=$(basename "$protein_dir"), Mutacao=$mutation_name, Ligante=$ligand_name"

            # Cria uma área de trabalho única
            WORK_DIR="$TEMP_BASE/$(basename "$protein_dir")_${mutation_name}_${ligand_name}"
            mkdir -p "$WORK_DIR/processed"

            # ----------------- PRÉ‑PROCESSAMENTO (dataset_creation) -----------------
            cd "$PPDOCK_REPO/dataset_creation"
            # Limpa resquícios de execuções anteriores
            rm -rf raw processed

            # Cria as pastas esperadas pelos scripts
            mkdir -p raw/test_data
            mkdir -p processed

            # Prepara os arquivos de entrada (nomes fixos exigidos pelos scripts)
            echo "$(basename "$protein_dir")_${mutation_name}_${ligand_name}" > raw/list.txt
            cp "$pdb_file"   "raw/test_data/$(basename "$protein_dir")_${mutation_name}.pdb"
            cp "$ligand_sdf" "raw/test_data/$(basename "$protein_dir")_${mutation_name}_${ligand_name}.sdf"

            # 1) Reordenar ligante (assumindo que o notebook foi convertido para .py)
            python 1_ligand_reorder.py || { echo "Erro no 1_ligand_reorder"; rm -rf raw processed; continue 2; }

            # 2) Processar proteína
            python 2_protein_process.py || { echo "Erro no 2_protein_process"; rm -rf raw processed; continue 2; }

            # 3) Processar ligante reordenado
            python 3_ligand_process.py || { echo "Erro no 3_ligand_process"; rm -rf raw processed; continue 2; }

            # 4) Criar .pt final
            python 4_data_pt.py || { echo "Erro no 4_data_pt"; rm -rf raw processed; continue 2; }

            # Copia o conteúdo gerado em processed para o WORK_DIR
            cp -r processed/* "$WORK_DIR/processed/"

            # Limpa os diretórios temporários no repositório
            rm -rf raw processed

            # ----------------- PREDIÇÃO DO BOLSÃO -----------------
            cd "$PPDOCK_REPO/1_Predict_Pocket"
            mkdir -p result
            ln -sfn "$WORK_DIR" dataset

            accelerate launch main_test.py \
                --data_path "$WORK_DIR" \
                --num_processes 1 \
            || { echo "Erro na predição do bolsão"; continue; }

            # Copia anchor.pt para processed (necessário para o docking)
            cp anchor.pt "$WORK_DIR/processed/"

            # ----------------- DOCKING -----------------
            cd "$PPDOCK_REPO/2_Dock"
            mkdir -p result
            ln -sfn "$WORK_DIR" dataset

            accelerate launch main_test.py \
                --data_path "$WORK_DIR" \
                --num_processes 1 \
            || { echo "Erro no docking"; continue; }

            # Coleta o resultado (ATENÇÃO: verifique o nome exato do arquivo de saída no código do main_test.py)
            if [ -f "docked_ligand.sdf" ]; then
                cp docked_ligand.sdf "$RESULT_DIR/$(basename "$protein_dir")_${mutation_name}_${ligand_name}.sdf"
                echo "OK: Docking salvo em $RESULT_DIR/$(basename "$protein_dir")_${mutation_name}_${ligand_name}.sdf"
            else
                echo "Aviso: arquivo de saída 'docked_ligand.sdf' não encontrado. Verifique o nome gerado pelo script."
            fi

            # Limpa a área de trabalho (comente para depurar)
            # rm -rf "$WORK_DIR"
        done
    done
done

conda deactivate
echo "Processamento concluído."
