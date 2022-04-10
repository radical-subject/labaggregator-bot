# Deploy

Есть виртуальный сервер с доступом по ssh

При коммите в ветку deploy GitHub Actions 
запускают скрипты которые обновляют git, 
копируют файлы на VPS и перезапукают докер


## Настройка виртуального сервера VPS 

### Ставим на сервер docker-compose

https://docs.docker.com/compose/install/

```
sudo curl -L "https://github.com/docker/compose/releases/download/1.29.2/docker-compose-$(uname -s)-$(uname -m)" -o /usr/local/bin/docker-compose
sudo chmod +x /usr/local/bin/docker-compose
docker-compose --version
```

https://docs.docker.com/engine/install/linux-postinstall/ 

```
sudo groupadd docker
sudo usermod -aG docker $USER
sudo systemctl enable docker.service
sudo systemctl enable containerd.service
```

## Настройка GitHub Actions

### Настраиваем Secrets

Для работы потребуется:
HOST - указываем ip адрес VPS
USERNAME - указываем имя пользователя VPS
SSHKEY - указываем приватный ssh ключ. его создание описано ниже

DEVELOPMENT_BOT_TOKEN - тестовый бот: используется для проверки работы всего на сервере
PRODUCTION_BOT_TOKEN - главный бот. публичный.

#### Делаем доступ по ssh на VPS

0. Заходим на VPS

`ssh <user>@<ip>`

1. Создаем ssh ключ для github

`ssh-keygen -t rsa -b 4096 -C "<email>"`

назовём файл github
без фразы without pass phrase

2. записываем публичный ключ в список доступа

`cat github.pub >> ~/.ssh/authorized_keys`

5. Печатаем и копируем в буфер приватный ключ

`cat github`

-----BEGIN OPENSSH PRIVATE KEY-----
b3BlbnNzaC1rZXktdjEAAAAABG5vbmUAAAAEbm9uZQAAAAAAAAABAAACFwAAAAdzc2gtcn
........
/gZe/RaaYlJdjQAAABRkb24tYW5kLWhvbWVAbWFpbC5ydQECAwQF
-----END OPENSSH PRIVATE KEY-----

6. Добавляем ключ в github secrets

`https://github.com/<username>/<repo>/settings/secrets/actions`

имя переменной: SSHKEY

