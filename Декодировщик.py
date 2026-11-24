# Декодер морзянки из probe.txt (FDTD-куб)
# Питон 3, нужна только numpy

import numpy as np

MORSE = {'.-':'A','-...':'B','-.-.':'C','-..':'D','.':'E','..-.':'F',
 '--.':'G','....':'H','..':'I','.---':'J','-.-':'K','.-..':'L','--':'M',
 '-.':'N','---':'O','.--.':'P','--.-':'Q','.-.':'R','...':'S','-':'T',
 '..-':'U','...-':'V','.--':'W','-..-':'X','-.--':'Y','--..':'Z',
 '-----':'0','.----':'1','..---':'2','...--':'3','....-':'4','.....':'5',
 '-....':'6','--...':'7','---..':'8','----.':'9','.-.-.-':'.'}

# ---------- 1. Загрузка ----------
d = np.loadtxt('probe.txt')
t, e = d[:,0], d[:,1]
t = t - t[0]                      # время от начала посылки

# ---------- 2. Огибающая: скользящее СКЗ ----------
# ВАЖНО: среднее НЕ вычитаем — морзянка записана в уровне поля,
# вычитание DC стирает сигнал (урок с первого прогона)
w = 100
env = np.sqrt(np.convolve(e**2, np.ones(w)/w, mode='same'))

# ---------- 3. Порог: посередине между плато посылки и нулём пауз ----------
thr = 0.4 * env.max()
on = env > thr

# ---------- 4. Серии одинаковых значений ----------
runs = []
cur, cnt = bool(on[0]), 1
for v in on[1:]:
    if v == cur: cnt += 1
    else: runs.append((cur, cnt)); cur, cnt = v, 1
runs.append((cur, cnt))

# отбрасываем краевые серии (до первой вспышки и после затухания)
first_on = next(i for i,(v,c) in enumerate(runs) if v)
last_on  = max(i for i,(v,c) in enumerate(runs) if v)
runs = runs[first_on:last_on+1]
if runs and not runs[0][0]: runs = runs[1:]
if runs and not runs[-1][0]: runs = runs[:-1]

# ---------- 5. Автокалибровка юнита ----------
# Ищем такой u, при котором все длины серий ближе всего к целым числам u.
# Это делает декодер устойчивым к искажениям губки и сглаживания:
# точки/паузы не обязаны быть ровно 1000 шагов.
lens = np.array([c for _,c in runs], float)
best_u, best_err = None, None
for u in np.linspace(0.5*lens.min(), 1.2*lens.min(), 400):
    r = lens/u
    err = np.sum((r - np.round(r))**2)
    if best_err is None or err < best_err:
        best_err, best_u = err, u
unit = best_u

# ---------- 6. Декодирование ----------
# Вспышка:  1 юнит -> '.', 2-3 юнита -> '-'
# Пауза:    1-2 юнита -> элемент внутри знака,
#           3-4 -> конец буквы, 5+ -> пробел (новое слово)
text, sym = '', ''
for v, c in runs:
    r = c/unit
    if v:
        sym += '.' if round(r) <= 1 else '-'
    else:
        rr = round(r)
        if rr <= 2:               # знак продолжается
            continue
        text += MORSE.get(sym, '?')
        sym = ''
        if rr >= 5:
            text += ' '
if sym:
    text += MORSE.get(sym, '?')

print('юнит =', round(unit), 'шагов')
print('РАСШИФРОВКА:', text.strip())