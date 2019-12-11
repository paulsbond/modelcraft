class Prune(Job):
    def __init__(self, cycle, xyzin, hklin, chains_only=False):
        super().__init__(cycle, "Coot pruning")
        self.xyzout = self.path("xyzout.pdb")
        self.run(xyzin, hklin, chains_only)
        self.finish()

    def run(self, xyzin, hklin, chains_only):
        coot_prune_dir = os.path.dirname(inspect.getfile(inspect.currentframe()))
        coot_prune_path = os.path.join(coot_prune_dir, "coot_prune.py")
        with open(coot_prune_path) as coot_prune:
            script = coot_prune.read()
        script += "deleted = prune(0, 1, 2, %s)\n" % chains_only
        script += "write_pdb_file(0, '%s')\n" % self.xyzout
        script += "import json\n"
        script += "json.dump(deleted, open('%s', 'w'), indent=2)\n" % self.path("deleted.json")
        script += "exit()"
        with open(self.path("script.py"), "w") as script_file:
            script_file.write(script)
        execute("coot", [xyzin, hklin,
                "--no-graphics", "--no-guano", "--script", self.path("script.py")])
        shutil.rmtree("coot-backup", ignore_errors=True)
        shutil.rmtree("coot-download", ignore_errors=True)
