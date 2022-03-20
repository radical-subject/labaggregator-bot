import clr
import os, sys
sys.path.append('../../../')
import modules.dbmodel as dbmodel
from System import DateTime
import dateparser
import pymorphy2

morph = pymorphy2.MorphAnalyzer()

dll_dir = "./Hors.0.9.30/lib/netstandard2.0/"
dllname = 'Hors'
path = r'%s%s' % (dll_dir, dllname)
sys.path.append(os.getcwd())
clr.AddReference(path)

from Hors import HorsTextParser
my_instance = HorsTextParser()

# input_text = "в прошлый четверг, Вождение машинки"
# user_id = 336091411

def hors_parse(user_input_text, user_id):
    result = my_instance.Parse(user_input_text, DateTime.Now)
    text = result.Text #-> string: "будет красивый закат"
    formatText = result.TextWithTokens #-> string: "{0} будет красивый закат"
    try:
        datefrom = result.Dates[0].DateFrom #-> DateTime: 26.03.2019 23:00:00
        dateto = result.Dates[0].DateTo #-> DateTime: 26.03.2019 23:00:00
    except IndexError:
        return None


    s = user_input_text.split()
    s_subtraction = text.split()

    list_of_categories = []
    for item in dbmodel.list_categories_including_archived(user_id):
        list_of_categories.append(item.lower())
    category = None
    for word in s_subtraction:
        if "," in list(word) or "." in list(word) or "!" in list(word):
            word = list(word)
            del word[-1]
            word = "".join(word)
        if morph.parse(word)[0].normal_form in list_of_categories:
            print ("категория определена автоматически")
            for i in range(len(list_of_categories)):
                if morph.parse(word)[0].normal_form == dbmodel.list_categories_including_archived(user_id)[i].lower():
                    category = dbmodel.list_categories_including_archived(user_id)[i]

    result_s = []
    for element in s:
        if element not in s_subtraction:
            result_s.append(element)
    s = " ".join(result_s)

    try:
        date = dateparser.parse(u'{}'.format(s), date_formats=['%d.%m.%Y', '%d-%m-%Y', '%d %m %Y', '%d/%m/%Y']).strftime("%Y-%m-%d %H:%M:%S")
        return [text, date, category]
    except AttributeError:
        if datefrom == None:
            return None
        else:
            return [text, dateparser.parse(u'{}'.format(str(datefrom))).strftime("%Y-%m-%d %H:%M:%S"), category]


# hors_parse(input_text)
