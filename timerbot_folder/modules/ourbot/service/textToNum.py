# %pip install pymorphy2
import pymorphy2
import re

morph = pymorphy2.MorphAnalyzer()

nums = {
    'полторы' : 1.5,
    'полтора' : 1.5,
    'четверть' : 0.25,
    'ноль' : 0,
    'одна' : 1,
    'две' : 2,
    'один' : 1,
    'два' : 2,
    'три' : 3,
    'четыре' : 4,
    'пять' : 5,
    'шесть' : 6,
    'семь' : 7,
    'восемь' : 8,
    'девять' : 9,
    'десять' : 10,
    'одиннадцать' : 11,
    'двенадцать' : 12,
    'тринадцать' : 13,
    'четырнадцать' : 14,
    'пятнадцать' : 15,
    'шестнадцать' : 16,
    'семнадцать' : 17,
    'восемнадцать' : 18,
    'девятнадцать' : 19,
    'двадцать' : 20,
    'тридцать' : 30,
    'сорок' : 40,
    'пятьдесят' : 50,
    'шестьдесят' : 60,
    'семьдесят' : 70,
    'восемьдесят' : 80,
    'девяносто' : 90,
    'сто' : 100,
    'двести' : 200,
    'триста' : 300,
    'четыреста' : 400,
    'пятьсот' : 500,
    'шестьсот' : 600,
    'семьсот' : 700,
    'восемьсот' : 800,
    'девятьсот' : 900,
}

THOUSANDS = {
    'тысяча' : 1000,
    'тысячи' : 1000,
    'тысяч' : 1000,
    2: ('миллион', 'миллиона', 'миллионов'),  # 10^6
    3: ('миллиард', 'миллиарда', 'миллиардов'),  # 10^9
    4: ('триллион', 'триллиона', 'триллионов'),  # 10^12
    5: ('квадриллион', 'квадриллиона', 'квадриллионов'),  # 10^15
    6: ('квинтиллион', 'квинтиллиона', 'квинтиллионов'),  # 10^18
    7: ('секстиллион', 'секстиллиона', 'секстиллионов'),  # 10^21
    8: ('септиллион', 'септиллиона', 'септиллионов'),  # 10^24
    9: ('октиллион', 'октиллиона', 'октиллионов'),  # 10^27
    10: ('нониллион', 'нониллиона', 'нониллионов'),  # 10^30
}


def str2num(text):
    words = text.split()
    result = 0
    result_list = []
    result_string = ''

    indexes_of_non_numerical_words = []
    for i in range(len(words)):
        if morph.parse(words[i])[0].tag.POS != "NUMR" and re.match(r'^([\s\d]+)$', words[i]) == None:
            if words[i] not in THOUSANDS.keys() and words[i] not in nums.keys():
                indexes_of_non_numerical_words.append(i)

    listik = []
    list_of_lists = []

    for i in range(len(words)):
        if i not in indexes_of_non_numerical_words:
            listik.append(words[i])
            if i == len(words)-1:
                list_of_lists.append(listik)
        else:
            if listik != []:
                list_of_lists.append(listik)
                list_of_lists.append(words[i])
                listik = []
            else:
                list_of_lists.append(words[i])

    final_list = []
    for element in list_of_lists:
        if str(type(element)) == "<class 'list'>":
            for word in element:
                if morph.parse(word)[0].tag.POS == "NUMR":
                    if word.lower() in nums.keys():
                        result += nums[word.lower()]
                elif word.lower() in nums.keys():
                    result += nums[word.lower()]
                elif word.lower() in THOUSANDS.keys():
                    if result == 0:
                        result = 1000
                    else:
                        result *= 1000

                elif re.match(r'^([\s\d]+)$', word) != None:
                    try:
                        result += float(word)
                    except ValueError:
                        if result == 0:
                            result_string += "1 " + word + " "
                        else:
                            result_string += str(result) + " " + word + " "
                            result = 0
            final_list.append(float(result))
            result = 0
        else:
            if morph.parse(element)[0].tag.POS != "NUMR" and element == list_of_lists[0]:
                if re.match(r's*(?:h|hour?|минута?|ч?|ч.?|час?|день?|секунда?|неделя?|неделю?)$', element) != None:
                    final_list.insert(0, 1)

    final_final_list = []
    i = 0
    for element in list_of_lists:
        if str(type(element)) == "<class 'list'>":
            final_final_list.append(final_list[i])
            i+=1
        else:
            if morph.parse(element)[0].tag.POS != "NUMR" and element == list_of_lists[0]:
                if re.match(r's*(?:h|hour?|минута?|ч?|ч.?|час?|день?|секунда?|неделя?|неделю?)$', element) != None:
                    final_final_list.append(final_list[i])
                    final_final_list.append(element)
                    i+=1
                else:
                    final_final_list.append(element)
            else:
                final_final_list.append(element)
    separator = ' '
    response_text = separator.join((str(i) for i in final_final_list))
    return response_text
