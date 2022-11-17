#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>

using namespace std;

// void read_record()
// {

// 	// File pointer
// 	fstream fin;

// 	// Open an existing file
// 	fin.open("test.csv", ios::in);

// 	// Get the roll number
// 	// of which the data is required
// 	long id, ids, count = 0;
// 	cout << "Enter the id number to display details: ";
// 	cin >> id;

// 	// Read the Data from the file
// 	// as String Vector
// 	vector<string> row;
// 	string line, word, temp;

// 	while (fin >> temp) {

// 		row.clear();

// 		// read an entire row and
// 		// store it in a string variable 'line'
// 		getline(fin, line);

// 		// used for breaking words
// 		stringstream s(line);

// 		// read every column data of a row and
// 		// store it in a string variable, 'word'
// 		while (getline(str, word, ', ')) {

// 			// add all the column data
// 			// of a row to a vector
// 			row.push_back(word);
// 		}

// 		// convert string to integer for comparision
// 		ids = stoi(row[0]);

// 		// Compare the roll number
// 		if (id == ids) {

// 			// Print the found data
// 			count = 1;
// 			cout << "Details of Roll " << row[0] << " : \n";
// 			cout << "Name: " << row[1] << "\n";
// 			cout << "Maths: " << row[2] << "\n";
// 			break;
// 		}
// 	}
// 	if (count == 0)
// 		cout << "Record not found\n";
// }

int main()
{
    string fname = "test.csv";

    vector<vector<string>> content;
    vector<string> row;
    string line, word;

    fstream file(fname, ios::in);
    if (file.is_open())
    {
        while (getline(file, line))
        {
            row.clear();

            stringstream str(line);

            while (getline(str, word, ','))
                row.push_back(word);
                cout << row[0] << endl;
                cout << row[1] << endl;
                cout << row[2] << endl;
            content.push_back(row);
        }
    }
    else {
        cout << "Could not open the file\n";
    }

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < content[i].size(); j++)
        {
            cout << content[i][j] << " ";
            
        }
        cout << "\n";
    }

    return 0;
}