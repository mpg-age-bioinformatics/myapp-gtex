# myapp-eg-container

Demo app for adding Apps as single containers to a myapp deployment.

docker-compose up -d
docker-compose exec server flask run --host 0.0.0.0 --port 8000
docker-compose exec gtex flask run --host 0.0.0.0 --port 8000

docker-compose exec gtex python3 -m smtpd -n -c DebuggingServer localhost:8025


