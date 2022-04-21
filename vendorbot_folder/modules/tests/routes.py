import logging
import json
import threading
import waitress
import werkzeug
from typing import List, Dict
from flask import Flask, request, make_response, send_file
from flask import json
from .core import core
from .storage import storage, FileBase
from .user import virtual_bot

logger = logging.getLogger(__name__)

app = Flask(__name__)

HOST = '127.0.0.1'
PORT = 5555
TELEGRAM_URL = f'http://{HOST}:{PORT}/'


@app.errorhandler(werkzeug.exceptions.BadRequest)
def handle_bad_request(e):
    logger.error('Server error: ' + str(e))
    return 'bad request!', 500


app.register_error_handler(500, handle_bad_request)


def result_ok(result: List, ok: bool = True) -> Dict:
    return {
            "ok": ok,
            "result": result
    }


def server_response(data):

    r = make_response(
       json.dumps(result_ok(data)),
       200
    )
    r.headers['Content-Type'] = 'application/json'
    return r


@app.route('/<token>/getMe', methods=['POST'])
def getMe(token: str):
    j = request.get_json()  # cause flask bug

    bot_id = int(token.split(':')[0])

    core.init_queue(bot_id)

    ret = virtual_bot.to_dict()
    return server_response(ret)


@app.route('/<token>/deleteWebhook', methods=['POST'])
def deleteWebhook(token: str):
    j = request.get_json()  # cause flask bug
    return server_response(True)


@app.route('/<token>/getUpdates', methods=['POST'])
def getUpdates(token: str):
    j = request.get_json()  # cause flask bug

    bot_id = int(token.split(':')[0])

    ret = core.get_updates(bot_id)

    return server_response(ret)


@app.route('/<token>/sendMessage', methods=['POST'])
def sendMessage(token: str):

    #bot_id = int(token.split(':')[0])

    data = request.get_json()

    chat_id = int(data['chat_id'])
    text = data['text']

    chat = {'id': chat_id, 'type': 'private'}  # simple chat structure

    ret = core.bot_send(bot_from=virtual_bot.to_dict(),
                        chat=chat,
                        text=text
                        )

    return server_response(ret)


@app.route('/<token>/getFile', methods=['POST'])
def getFile(token: str):

    data = request.get_json()
    file_id = data['file_id']

    document = storage.get(file_id)
    if document:
        file = FileBase(document,
                        file_path=f'file/{file_id}')

        return server_response(file.to_dict())

    return result_ok([], False)


@app.route('/<token>/file/<file_id>', methods=['GET'])
def file(token: str, file_id: str):

    document = storage.get(file_id)

    return send_file(document.local_path, as_attachment=True)


def start_server():
    s = waitress.create_server(app, host=HOST, port=PORT)

    def run(s):
        try:
            s.run()
        except Exception as err:
            logger.error(err)
        del s
        logger.info("server deleted")

    t = threading.Thread(target=run, args=(s,))
    t.daemon = True
    t.start()
    return s, t


def shutdown_server(s, t):
    logger.info('start close server')
    s.close()
    logger.info('server closed')

    if t.is_alive():
        logger.info('start join thread')
        t.join(timeout=2.0)
        logger.info('thread joined')


if __name__ == "__main__":

    from time import sleep

    s, t = start_server()

    while True:
        sleep(1)

    shutdown_server(s, t)
