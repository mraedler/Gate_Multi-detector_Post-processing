#include "../include/utils.h"

#include <iostream>
#include <TFile.h>
#include <TList.h>
#include <TKey.h>
#include <TObjArray.h>
#include <TSystemDirectory.h>
#include <TObjString.h>
#include <sys/wait.h>
#include <unistd.h>

void listAvailableTrees(const TFile* file) {
    TList* keys = file->GetListOfKeys();
    TIter next(keys);
    TKey* key;
    while ((key = (TKey*)next())) {
        TObject* obj = key->ReadObj();
        if (obj->InheritsFrom(TTree::Class())) {
            std::cout << "Found tree: " << obj->GetName() << std::endl;
        }
    }
    std::cout << std::endl;
}



void listAvailableBranches(TTree* tree) {
    TObjArray* branches = tree->GetListOfBranches();
    for (int ii = 0; ii < branches->GetEntries(); ii++) {
        TBranch* branch = (TBranch*)branches->At(ii);
        std::cout << "Found branch: " << branch->GetName() << " (title: " << branch->GetTitle() << ")" << std::endl;
    }
    std::cout << std::endl;
}



void listAvailableLeaves(TTree* tree) {
    TObjArray* leaves = tree->GetListOfLeaves();
    for (int ii = 0; ii < leaves->GetEntries(); ii++) {
        TLeaf* leaf = (TLeaf*)leaves->At(ii);
        std::cout << "Found leaf: " << leaf->GetName() << " (type: " << leaf->GetTypeName() << ")" << std::endl;
    }
    std::cout << std::endl;
}



bool checkIfLeafExists(TTree* tree, const std::string& checkLeafName) {
    bool leafExists = false;
    TObjArray* leaves = tree->GetListOfLeaves();
    for (int ii = 0; ii < leaves->GetEntries(); ++ii) {
        std::string leafName(leaves->At(ii)->GetName());

        if (leafName.compare(0, checkLeafName.length(), checkLeafName) == 0) {
            //std::cout << "Found leaf: " << leafName << std::endl;
            leafExists = true;
        }
    }

    return leafExists;
}



TFile* getTFile(const TString& fullPath, const char* mode, const bool verbose) {
    TFile* file = TFile::Open(fullPath, mode);

    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file " << fullPath << std::endl;
        std::exit(1);
    }
    if (verbose) {listAvailableTrees(file);}

    return file;
}



TTree* getTTree(TFile* file, const TString& treeName, const bool verbose) {
    auto* tree = file->Get<TTree>(treeName);
    if (!tree) {
        std::cerr << "Tree '" << treeName << "' not found." << std::endl;
        std::exit(1);
    }
    //if (verbose) {listAvailableBranches(tree);}
    if (verbose) {listAvailableLeaves(tree);}

    return tree;
}



void runInSeparateProcesses(const size_t nJobs, const JobFn job, const int maxProcesses) {
    std::vector<pid_t> children;

    for (size_t i = 0; i < nJobs; ++i) {
        pid_t pid = fork();

        if (pid == 0) {
            job(i);
            _exit(0);
        }

        if (pid > 0) {
            children.push_back(pid);
            if ((int)children.size() >= maxProcesses) {
                wait(nullptr);
                children.pop_back();
            }
            continue;
        }

        perror("fork");
        std::exit(1);
    }

    while (wait(nullptr) > 0);
}



void runSequentially(const size_t nJobs, const JobFn job) {
    for (unsigned int ii = 0; ii < nJobs; ii++) {
        job(ii);
    }
}



std::vector<TString> getListOfRootFilePaths(const char* dirPath, const bool verbose) {
    // Get list of .root files
    const TSystemDirectory dir("dir", dirPath);
    const TList* files = dir.GetListOfFiles();

    // Allocate the array of .root files
    std::vector<TString> fullPaths;

    if (!files) {
        std::cerr << "Error: cannot open directory " << dirPath << std::endl;
        return fullPaths;
    }

    // files->Sort();

    TString dirPathTString(dirPath);
    if (!dirPathTString.EndsWith("/")) {
        dirPathTString += "/";
    }

    TIter next(files);
    TSystemFile* file;
    while ((file = dynamic_cast<TSystemFile*>(next()))) {
        TString fileName = file->GetName();
        if (!file->IsDirectory() && fileName.EndsWith(".root")) {
            fullPaths.push_back(dirPathTString + fileName);
        }
    }

    if (verbose) {
        for (const auto & fullPath : fullPaths) {
            std::cout << fullPath << std::endl;
        }
    }

    return fullPaths;
}



std::vector<TString> splitPath(const TString& path) {
    TObjArray* parts = path.Tokenize("/");

    // directories (excluding filename)
    std::vector<TString> dirs;
    for (int ii = 0; ii < parts->GetEntries(); ++ii) {
        dirs.emplace_back(((TObjString*)parts->At(ii))->GetString());
    }

    return dirs;
}



TString assignGantryName(Int_t gantryID1, Int_t gantryID2) {
    bool isBI1 = (gantryID1 == 2);
    bool isBI2 = (gantryID2 == 2);

    TString gantryName;
    if (!isBI1 && !isBI2) {
        gantryName = "TB-TB";
    } else if (isBI1 ^ isBI2) {
        gantryName = "TB-BI";
    } else if (isBI1 && isBI2) {
        gantryName = "BI-BI";
    } else {
        gantryName = "Unknown";
        std::cerr << "Warning: unknown gantry.\n";
    }
    return gantryName;
}
