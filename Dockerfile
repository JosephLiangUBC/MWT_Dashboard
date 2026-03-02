# app/Dockerfile

FROM python:3.11-slim

WORKDIR /MWT_Dashboard

COPY requirements.txt ./requirements.txt

RUN pip3 install -r requirements.txt

EXPOSE 8502

COPY . /MWT_Dashboard

ENTRYPOINT ["streamlit", "run", "MWT_dashboard.py", "--server.port=8502"]