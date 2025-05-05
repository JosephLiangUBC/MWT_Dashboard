# MWT_Dashboard
 Data Dashboard for MWT

Youtube Overview:

[![MWT Dashboard Introduction Video](http://img.youtube.com/vi/xbQ0CRUCnZs/0.jpg)](https://youtu.be/xbQ0CRUCnZs "MWT Dashboard Introduction Video")



### Installing Dependencies

#### Using pip:

1. Clone the repository
2. Create and activate a vitual environemnt
```
python -m venv rankinlab
source rankinlab/bin/activate     # macOS/Linux
rankinlab\Scripts\activate.bat    # Windows
```

3. Install the required dependencies
```
pip install -r requirements.txt

```
4. Update `requirements.txt` when new packages are installed 
```
pip freeze > requirements.txt`
```

#### Using conda 

1. Clone the repository
2. Create and activate conda environment
```
conda env create -f environment.yml
conda activate rankinlab
```
3. Update `environment.yml` when new packages are installed
```
conda env update --file environment.yml --prune
```