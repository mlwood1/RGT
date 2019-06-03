''' docstring'''
import xlwt

class ExcelWriter():
    """docstring for ClassName"""
    def __init__(self):
        self.wb = xlwt.Workbook()


    def add_table_to_sheet(self, table, sheet_name):
        row = 0
        ws = self.wb.add_sheet(sheet_name)
        for key in table.keys():
            ws.write(row,0,key)
            values = table[key]
            if isinstance(values, list): #check if the value is a single value or a list
                for idx, val in enumerate(values):
                    print(idx, val)
                    ws.write(row,1+idx, val)
            else:
                    ws.write(row,1, table[key])

            row+=1
    
    def save_file(self, file_name="DefaultName.xls"):
        self.wb.save(file_name)

    @staticmethod	
    def write_to_excel(table, sheet_name="DefaultName", file_name="DefaultName.xls"):
        row = 0
        wb = xlwt.Workbook() #Creat a new workbook
        ws = wb.add_sheet(sheet_name)
        for key in table.keys():
            ws.write(row,0,key)
            ws.write(row,1, table[key])
            row+=1
        wb.save(file_name)