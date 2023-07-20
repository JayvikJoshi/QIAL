
import os
import openpyxl

def collect_filenames_to_files(directory_path, text_file, excel_file):
    try:
        # Get the list of filenames in the directory
        filenames = os.listdir(directory_path)

        # Write filenames to the text file
        with open(text_file, 'w') as text_file:
            for filename in filenames:
                text_file.write(filename + '\n')

        # Write filenames to the Excel file
        workbook = openpyxl.Workbook()
        sheet = workbook.active
        for index, filename in enumerate(filenames, start=1):
            sheet.cell(row=index, column=1, value=filename)

        workbook.save(excel_file)

        print("Filenames collected and written to", text_file, "and", excel_file)
    except Exception as e:
        print("An error occurred:", str(e))

# Replace 'directory_path' with the path of the directory you want to collect filenames from
# Replace 'output_text_file.txt' with the desired name of the output text file
# Replace 'output_excel_file.xlsx' with the desired name of the output Excel file
collect_filenames_to_files('/Users/Jayvik/Desktop/All_Data/processed_t1', '/Users/Jayvik/Desktop/filenames_txt.txt', '/Users/Jayvik/Desktop/filenames_excel.xlsx')
