#include <cassert>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <io/path_list.hpp>

namespace fs = std::filesystem;

int main()
{
    const fs::path dir =
        fs::temp_directory_path() / "agnsf_test_path_list";

    fs::remove_all(dir);
    fs::create_directories(dir);

    const fs::path file = dir / "paths.txt";

    {
        std::ofstream out(file);
        assert(out);

        out << "# light-curve files\n";
        out << "\n";
        out << "  /data/lc1.csv  \n";
        out << "/data/lc2.fits\n";
        out << "# a comment\n";
        out << "/data/lc3.fits.gz\n";
    }

    const std::vector<std::string> paths =
        agnsf::io::read_path_list(file.string());

    assert(paths.size() == 3);
    assert(paths[0] == "/data/lc1.csv");
    assert(paths[1] == "/data/lc2.fits");
    assert(paths[2] == "/data/lc3.fits.gz");

    bool thrown = false;

    try {
        agnsf::io::read_path_list(
            (dir / "missing.txt").string()
        );
    }
    catch (const std::runtime_error&) {
        thrown = true;
    }

    assert(thrown);

    fs::remove_all(dir);

    return 0;
}
