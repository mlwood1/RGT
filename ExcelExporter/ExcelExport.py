''' docstring'''
from openpyxl import Workbook

class ExcelWriter():
    """docstring for ClassName"""
    def __init__(self):
        #self.wb = xlwt.Workbook()
        self.wb = Workbook()
        std = self.wb.get_sheet_by_name('Sheet')
        self.wb.remove_sheet(std)

    def add_table_to_sheet(self, table, sheet_name):
        row = 1
        ws = self.wb.create_sheet(sheet_name)
        ws.title = sheet_name

        for key in table.keys():
            ws.cell(row=row, column=1, value=key)
            values = table[key]
            if isinstance(values, list): #check if the value is a single value or a list
                for idx, val in enumerate(values):
                    ws.cell(row=row, column=2+idx, value=val)

            else:   
                    ws.cell(row=row, column=2, value=table[key])

            row+=1
    
    def save_file(self, file_name="DefaultName.xlsx"):
        self.wb.save(file_name)

    @staticmethod	
    def write_to_excel(table, sheet_name="DefaultName", file_name="DefaultName.xls"):
        row = 1
        wb = Workbook()
        ws = wb.create_sheet(sheet_name)
        for key in table.keys():
            ws.cell(row=row, column=1, value=key)
            ws.cell(row=row, column=2, value=table[key])
            row+=1
        wb.save(file_name)