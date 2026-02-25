# app/Dockerfile

FROM python:3.11-slim

WORKDIR /MWT_Dashboard

COPY requirements.txt ./requirements.txt

RUN pip3 install -r requirements.txt

EXPOSE 8501

HEALTHCHECK CMD curl --fail http://localhost:8501/_stcore/health

COPY . /MWT_Dashboard

ENTRYPOINT ["streamlit", "run", "MWT_dashboard.py", "--server.port=8501"]