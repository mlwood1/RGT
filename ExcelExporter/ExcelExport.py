''' docstring'''
import xlwt

class ExcelWriter(object):
    """docstring for ClassName"""
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