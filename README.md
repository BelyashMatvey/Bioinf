# Модуль Bioinf
## Аннотация
Модуль **Bioinf** представляет содержит все необходимые для начинающего биоинформатика функции.
Он позволяет базово преобразовывать последовательности нуклеиновых кислот, а также фильтровать 
данные в формате _fastq_.

## Оглавление
- [Основные функции](#main-functions)
- [Используемые данные](#data)
- [Авторы](#authors)

<a name="main-functions"><h2>Основные функции</h2></a>
### filter_fastq
**filter_fastq** принимает на вход словарь с данными прочтения последовательности нуклеотидов и 
и параметры для их фильтрации (GC состав, длина прочтения, среднее качество прочтения). 

Функция возвращает словарь с отфильтрованными по соответствующим
параметрам данными.

### run_dna_rna_tools
**run_dna_rna_tools** принимает на вход несколько строк нуклеиновых кислот и операцию, 
которую необходимо произвести над ними. Возвращает результат операции для каждой заданной последовательности
Функция использует пакет _modules_, содержащий модули для соответствующих операций. Если операция не может быть
применена данной нуклеотидной последовательности или входная строка не является нуклеотидной последовательностью -
выводится соответствующая ошибка.
#### Возможные операции:
 - **complement** - возвращает комплементарную данной последовательность
 - **reverse** - возвращает развернутую последовательность
 - **transcribe** - возвращает транскрибированную строку
 - **reverse_complement** - возвращает обратную комплементарную последовательность
 - **transcribe_dna_complement** - возвращает комплементарную последовательность для транскрипции заданной последовательности

<a name="data"><h2>Используемые данные</h2></a>
В пакете _modules_ содержится модуль **const_dicts_and_strs**, который содержит в себе необходимые словари для операций над
последовательностями, а также тексты возможных ошибок.

<a name="authors"><h2>Авторы</h2></a>
### **Беляков Матвей** 

