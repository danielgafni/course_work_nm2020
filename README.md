# Установка

**Следующие команды выполняются из терминала в директории course\_work\_nm\_2020.**

При помощи [conda](https://www.anaconda.com/distribution/):

```shell
conda create -y --name cwnm2020
conda activate cwnm2020
conda install -y -c plotly -c numpy numba plotly plotly-orca jupyter tqdm psutil requests
```

# Воспроизведение

Находясь в той же директории выполните из терминала:

```shell
jupyter notebook
```

Откроется Jupyter Notebook. Выберите файл `course_work_nm_2020.ipynb`.

# Просмотр готового .ipynb при помощи веб-браузера

Работу можно просмотреть без необходимости установки библиотек и запуска Jupyter Notebook сервера. Для этого перейдите по ссылке

[https://nbviewer.jupyter.org/github/danielgafni/course_work_nm2020/blob/master/Heat_equation_2D_interactive.ipynb](https://nbviewer.jupyter.org/github/danielgafni/course_work_nm2020/blob/master/Heat_equation_2D_interactive.ipynb)

Также интерактивные графики доступны в виде `.html` файлов в папке `results` после первого запуска `Heat_equation_2D.ipynb`.
