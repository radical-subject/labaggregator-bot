import os
try: 
    f = open("./secrets/token.secret", "r")
    token = f.read()

except:
    token = os.environ["BOT_TOKEN"]