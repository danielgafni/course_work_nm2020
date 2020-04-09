# Установка

**Следующие команды выполняются из терминала в директории course\_work\_nm\_2020.**

При помощи [conda](https://www.anaconda.com/distribution/):

```shell
conda create -y --name cwnm2020

conda install -y plotly jupyter numba tqdm

conda activate cwnm2020
```

# Воспроизведение

Находясь в той же директории выполните из терминала:

```shell
jupyter notebook
```

Откроется Jupyter Notebook. Выберите файл `course_work_nm_2020.ipynb`.

# Просмотр готового .ipynb при помощи веб-браузера

Работу можно просмотреть без необходимости установки библиотек и запуска Jupyter Notebook сервера. Для этого перейдите по ссылке

[ссылка на nbviewer](asdadas)

Также графики доступны в виде .html файлов в папке `results`.

Файл `course_work.pdf` получен с помощью [doPDF](https://ru.dopdf.com/) из .ipynb файла.