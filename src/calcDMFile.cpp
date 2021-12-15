#include "dm_compute.h"
#include <chrono>
#include <thread>
#include <unistd.h>
#include <limits.h>

/*
Randomly sample mssm parameter space and calculate Higgs masses.

TO DO:
1. Fix the seeding issues, it is possible that the initial seed is the same due to being distributed acrossed different nodes.
Last Updated: 10/24/2017
*/

vector<string> readConfig()
{
	vector<string> parameters = {};
	ifstream configFile;
	// This .cpp file uses relative paths!
	configFile.open("../config.txt");
	if(!configFile)
	{
		cout << "The parent directory must contain a file config.txt." << endl;
		exit(1);
	}
	string line;
	while(getline(configFile, line))
		if(line.length()!= 0 && line[0] != '#')
		{
			string indicator = line.substr(0,2);
			vector<string> split_content;
			string content;
			boost::split(split_content, line, boost::is_any_of(" "))[1];
			content = split_content[1];
			boost::trim(content);
			if(indicator != "##")
                                parameters.push_back(content);
		}
	configFile.close();
	return parameters;
}

int main(int argc, char* argv[])
{
	// Read in parameters from config file
	vector<string> parameters = readConfig();
	string micrOmegasDir = parameters[0];
        string softSUSYDir = parameters[1];
	string param_space = parameters[2];
        string unranFilesDir = parameters[3];
        string outputFilesDir = parameters[4];
	string temporaryDir_root = parameters[5];
	int number_threads = stoi(parameters[6]);
        // Get the batch ID number to ensure unique seeds
        int fileNum;
        if(argc < 1)
	{
	    cout << "You need to pass a seed with this program" << endl;
            return 0;
	}
        else
            fileNum = atoi(argv[1]);
	// Read in points from file.
	stringstream ss;
    ifstream readfile(unranFilesDir + "/generated_points_" + to_string(fileNum) + ".txt");
    cout << unranFilesDir + "/generated_points_" + to_string(fileNum) + ".txt" << endl;
    string temporaryDir = temporaryDir_root + "/" + "working_directory_"+to_string(fileNum);
    string create_temporaryDir = string("mkdir " + temporaryDir);
    system(create_temporaryDir.c_str());
    DM dm(micrOmegasDir, softSUSYDir, temporaryDir);
    vector<vector<string>> cMSSM_points;
    string line;
    while(getline(readfile, line))
    {
        vector<string> one_point;
        boost::split(one_point, line, boost::is_any_of(","));
        cMSSM_points.push_back(one_point);
    } 
    int start_s=time(NULL);
    cout << "Beginning while loop" << endl;
    #pragma omp parallel for num_threads(number_threads)
    for(int i=0; i < cMSSM_points.size(); i++)
    {
        string thread_id = to_string(omp_get_thread_num());
	double sgnmu;
	double m0, m12, a0, tanb, m1, m2, m3, at, ab, ata, mel, mtal, mer, mtar, mq1l, mq3l, mur, mtr, mdr, mbr, mmu, ma;
        int clock = time(0);
	ofstream thread_file;	
 	thread_file.open(temporaryDir + "/outputjh" + thread_id + ".txt", ios_base::app);
	vector<double> mssm_point;
	string result; string fail_key;
        m0 = atof(cMSSM_points[i][0].c_str());
        m12 = atof(cMSSM_points[i][1].c_str());
        a0 = atof(cMSSM_points[i][2].c_str());
        tanb = atof(cMSSM_points[i][3].c_str());
        if(atof(cMSSM_points[i][4].c_str()) > 0)
            sgnmu = 1;
        else
            sgnmu = -1;
        // Value to return in case of failed theory (just call everything -1).
        // Running test.o will give a copy of the fail key that can be copy + pasted.
        fail_key = "-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1";
	if(param_space == "2D")
	{
	    mssm_point = {m0, m12, sgnmu};
	    result = dm.run_calc(mssm_point, thread_id);
            if(result != "Error")
		thread_file << m0 << "," << m12 << "," << result << endl;
            else
                thread_file << m0 << "," << m12 << "," << fail_key << endl;		
	}			
	if(param_space == "cMSSM")
	{
	    mssm_point = {m0, m12, a0, tanb, sgnmu};
	    result = dm.run_calc(mssm_point, thread_id);
            if(result != "Error")			
	        thread_file << m0 << "," << m12 << "," << a0 << "," << tanb << "," << sgnmu << "," << result << endl;
            else
                thread_file << m0 << "," << m12 << "," << a0 << "," << tanb << "," << sgnmu << "," << fail_key << endl;
	}
	if(param_space == "pMSSM")
	{   
            m1 = atof(cMSSM_points[i][0].c_str());
            m2 = atof(cMSSM_points[i][1].c_str());
            m3 = atof(cMSSM_points[i][2].c_str());
            mmu = atof(cMSSM_points[i][3].c_str());
            ma = atof(cMSSM_points[i][4].c_str());
            at = atof(cMSSM_points[i][5].c_str());
            ab = atof(cMSSM_points[i][6].c_str());
            ata = atof(cMSSM_points[i][7].c_str());
            mel = atof(cMSSM_points[i][8].c_str());
            mtal = atof(cMSSM_points[i][9].c_str());
            mer = atof(cMSSM_points[i][10].c_str());
            mtar = atof(cMSSM_points[i][11].c_str());
            mq1l = atof(cMSSM_points[i][12].c_str());
            mq3l = atof(cMSSM_points[i][13].c_str());
            mur = atof(cMSSM_points[i][14].c_str());
            mtr = atof(cMSSM_points[i][15].c_str());
            mdr = atof(cMSSM_points[i][16].c_str());
            mbr = atof(cMSSM_points[i][17].c_str());
            tanb = atof(cMSSM_points[i][18].c_str());
    	    mssm_point = {m1, m2, m3, mmu, ma, at, ab, ata, mel, mtal, mer, mtar, mq1l, mq3l, mur, mtr, mdr, mbr, tanb};
	    result = dm.run_calc(mssm_point, thread_id);
	    if(result != "Error")
	 	thread_file << m1 << "," << m2 << "," << m3 << "," << mmu << "," << ma << "," << at << "," << ab << "," << ata << "," << 
					mel << "," << mtal << "," << mer << "," << mtar << "," << mq1l << "," << mq3l << "," << mur << "," << mtr << "," <<
					mdr << "," << mbr << "," << tanb << "," << result << endl;
            else
                thread_file << m1 << "," << m2 << "," << m3 << "," << mmu << "," << ma << "," << at << "," << ab << "," << ata << "," << 
					mel << "," << mtal << "," << mer << "," << mtar << "," << mq1l << "," << mq3l << "," << mur << "," << mtr << "," <<
					mdr << "," << mbr << "," << tanb << "," << fail_key << endl;
	}		
    }
    // Combine all thread files into one node file.
    ofstream writefile;
    writefile.open(outputFilesDir + "/ran_points_" + to_string(fileNum)+ ".txt", ios_base::app);
    for(int i=0; i<number_threads; i++)
    {
	string line;
	ifstream readfile(temporaryDir + "/outputjh" + to_string(i) + ".txt");		
	while(getline(readfile, line))
		writefile << line << endl; 
	readfile.close();
	string cmd = string("rm " + temporaryDir + "/outputjh" + to_string(i) + ".txt"); 
	system(cmd.c_str());
    }	
    writefile.close();
    int stop_s=time(NULL);
    cout << "time per point: " << double(stop_s-start_s) / cMSSM_points.size() << " seconds" << endl;
}
