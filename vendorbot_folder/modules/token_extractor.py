import os

#TODO удалить - замена файл .env в корне с BOT_TOKEN=....
try: 
    f = open("./secrets/token.secret", "r")
    token = f.read()

except:
    token = os.environ["BOT_TOKEN"]