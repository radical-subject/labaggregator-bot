import sys
sys.path.append('..')
import modules.pytimeparse as pytimeparse
import modules.textToNum as textToNum

def parse_text_to_time_in_mins(user_input_text):

    try:
        minutes = pytimeparse.timeparse(textToNum.str2num(user_input_text))/60
        response = float("{:.1f}".format(minutes))
        return response
    except:
        return None
