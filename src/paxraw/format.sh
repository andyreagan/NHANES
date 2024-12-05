# BUILD_JUPYTERLAB_EXTENSION=1 pip install git+https://github.com/mwouts/jupytext.git

format () {
    black *.py --line-length=100
    isort *.py -l 100 -p utils --profile black  # or set config to ../../.isort.cfg
    flake8 *.py --config=../../.flake8
}

update_notebooks () {
    NOTEBOOKS=($(ls *.ipynb))

    format

    for FILE in $NOTEBOOKS; do
        BASE=$(basename $FILE .ipynb)
        jupytext $BASE.py --to ipynb
    done
}

update_scripts () {
    NOTEBOOKS=($(ls *.ipynb))

    for FILE in $NOTEBOOKS; do
        BASE=$(basename $FILE .ipynb)
        jupytext $BASE.ipynb --to py
    done

    format

    for FILE in $NOTEBOOKS; do
        BASE=$(basename $FILE .ipynb)
        jupytext $BASE.py --to ipynb --update
    done
}
