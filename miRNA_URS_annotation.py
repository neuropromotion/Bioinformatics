import requests  
from bs4 import BeautifulSoup
import re

def get_human_mirna_name(urs_code):
    """Извлекает название человеческой микроРНК по URS коду"""
    url = f'https://rnacentral.org/rna/{urs_code}/9606'
    
    response = requests.get(url)
    if response.status_code == 200:
        soup = BeautifulSoup(response.text, 'html.parser')
        text_content = soup.get_text()
        
        # Паттерн для поиска названия после "Homo sapiens (human)"
        pattern = r'Homo sapiens \(human\)\s+([^\n\r\|]+?)(?:\s+\||\n|\r|$)'
        matches = re.findall(pattern, text_content)
        
        if matches:
            human_mirna = matches[0].strip()
            # Извлекаем только само название (убираем описание в скобках)
            mirna_name = human_mirna.split('(')[0].strip()
            return mirna_name
        else:
            # Альтернативные паттерны для поиска
            broader_patterns = [
                r'(Hsa-Mir-[^\s\n\r\|]+)',
                r'(hsa-miR-[^\s\n\r\|]+)',
                r'Homo sapiens.*?([Hh]sa-[Mm]i[rR]-[^\s\n\r\|]+)'
            ]
            
            for pattern in broader_patterns:
                matches = re.findall(pattern, text_content, re.IGNORECASE)
                if matches:
                    return matches[0]
            
            return ''
    else:
        return 'Error'

# Пример использования
urs_code = 'URS000023B77E'
result = get_human_mirna_name(urs_code)
print(f"Результат для {urs_code}: {result}")
#Результат для URS000023B77E: hsa-miR-200a-5p
