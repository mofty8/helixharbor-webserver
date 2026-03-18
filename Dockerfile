FROM python:3.9.17-bullseye

COPY requirements.txt /

EXPOSE 5005

RUN pip3 install --upgrade pip \
 && pip3 install --default-timeout=300 --no-cache-dir -r /requirements.txt

COPY . /app

WORKDIR /app

ENTRYPOINT ["./gunicorn.sh"]



