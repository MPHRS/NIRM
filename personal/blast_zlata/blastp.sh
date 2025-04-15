#!/bin/bash

# Проверка, что скрипт запускается из правильной папки (data)
if [ ! -d "data" ]; then
  echo "Ошибка: Папка 'data' не найдена!"
  exit 1
fi

# Создаем папку для результатов, если её нет
mkdir -p blast_results

# Получаем список всех папок с образцами в папке "data"
organisms=($(ls -d data/*/))

# Создание BLAST-баз данных для каждого образца
echo "Создание BLAST баз данных для всех образцов..."
for org in "${organisms[@]}"; do
    # Проверка наличия файла protein.faa в каждой папке
    if [ -f "${org}protein.faa" ]; then
        echo "Создаем базу данных для ${org}..."
        makeblastdb -in "${org}protein.faa" -dbtype prot -out "${org}protein_db"
    else
        echo "❌ Ошибка: Файл ${org}protein.faa не найден!"
    fi
done

# Запуск BLASTP для попарного сравнения всех образцов
echo "Запуск BLASTP для попарного сравнения..."
for org1 in "${organisms[@]}"; do
    for org2 in "${organisms[@]}"; do
        if [ "$org1" != "$org2" ]; then
            # Убираем слэши из имен папок
            org1_name=$(basename "$org1")
            org2_name=$(basename "$org2")

            # Проверяем, существуют ли файлы protein.faa для обоих организмов
            if [ -f "${org1}protein.faa" ] && [ -f "${org2}protein.faa" ]; then
                echo "Сравниваем $org1_name с $org2_name..."
                blastp -query "${org1}protein.faa" -db "${org2}protein_db" -out "blast_results/${org1_name}_${org2_name}.txt" -evalue 1e-5 -outfmt 6 -num_threads 4
            fi
        fi
    done
done

echo "✅ BLASTP анализ завершен. Результаты в папке blast_results."

