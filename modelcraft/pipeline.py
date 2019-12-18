from modelcraft.arguments import parse
from modelcraft.buccaneer import Buccaneer
from modelcraft.findwaters import FindWaters
from modelcraft.gemmineer import model_stats
from modelcraft.prune import Prune
from modelcraft.refmac import Refmac
import shutil


class Pipeline():
    def __init__(self, argument_list):
        print("# ModelCraft")
        print("\nPlease cite [paper to be published]")
        self._args = parse(argument_list)
        self._initialise()
        self._run()

    def _initialise(self):
        self.min_rwork = 1
        self.min_rfree = 1
        self.min_fragments_built = 999
        self.max_longest_fragment = 0
        self.max_residues_built = 0
        self.max_residues_sequenced = 0
        self.cycles_without_improvement = 0
        self.jobs = {}

    def _run(self):
        args = self._args
        for self.cycle in range(1, args.cycles + 1):
            print("\n## Cycle %d\n" % self.cycle)
            self.jobs[self.cycle] = []
            if self.cycle == 1:
                if args.colin_hl is None and args.colin_phifom is None and args.mr_model is not None:
                    refmac = self._refmac(args.mr_model)
                    buccaneer = self._buccaneer(refmac.hklout, args.xyzin)
                else:
                    buccaneer = self._buccaneer(args.hklin, args.xyzin)
            else:
                coot = self._prune(refmac.xyzout, refmac.hklout)
                refmac = self._refmac(coot.xyzout)
                buccaneer = self._buccaneer(refmac.hklout, refmac.xyzout)
            refmac = self._refmac(buccaneer.xyzout)
            coot = self._prune(refmac.xyzout, refmac.hklout, chains_only=True)
            if args.add_waters and self.min_rwork < 0.4:
                coot = self._findwaters(coot.xyzout, refmac.hklout)
            refmac = self._refmac(coot.xyzout)
            self.process_cycle_output(refmac)
            if self.cycle > 1:
                self.remove_job_directories(self.cycle - 1)
            if args.auto_stop and self.cycles_without_improvement == 4:
                break
        self.remove_job_directories(self.cycle)

    def _job_directory(self, name):
        return "%02d.%02d_%s" % (self.cycle, len(self.jobs[self.cycle]) + 1, name)

    def _buccaneer(self, hklin, xyzin):
        directory = self._job_directory("buccaneer")
        cycles = 3 if self.cycle == 1 else 2
        job = Buccaneer(self._args, directory, hklin, xyzin, cycles)
        self.jobs[self.cycle].append(job)
        return job

    def _refmac(self, xyzin):
        directory = self._job_directory("refmac")
        use_phases = self._args.unbiased and self.min_rwork > 0.35
        job = Refmac(self._args, directory, xyzin, use_phases)
        self.jobs[self.cycle].append(job)
        return job

    def _prune(self, xyzin, hklin, chains_only=False):
        directory = self._job_directory("prune_chains" if chains_only else "prune")
        job = Prune(directory, xyzin, hklin, chains_only)
        self.jobs[self.cycle].append(job)
        return job

    def _findwaters(self, xyzin, hklin):
        directory = self._job_directory("findwaters")
        job = FindWaters(directory, xyzin, hklin)
        self.jobs[self.cycle].append(job)
        return job

    def improved(self, cycle):
        required_improvement = 0.02
        improvement = (self.min_rwork - cycle["rwork"]) / self.min_rwork
        if improvement > required_improvement:
            return True
        improvement = (cycle["residues_built"] - self.max_residues_built) / float(self.max_residues_built)
        if improvement > required_improvement:
            return True
        improvement = (cycle["residues_sequenced"] - self.max_residues_sequenced) / float(self.max_residues_sequenced)
        if improvement > required_improvement:
            return True
        improvement = (self.min_fragments_built - cycle["fragments_built"]) / float(self.min_fragments_built)
        if improvement > required_improvement:
            return True
        improvement = (cycle["longest_fragment"] - self.max_longest_fragment) / float(self.max_longest_fragment)
        if improvement > required_improvement:
            return True
        return False

    def process_cycle_output(self, refmac):
        cycle = {
            "cycle": self.cycle,
            "rwork": refmac.final_rwork,
            "rfree": refmac.final_rfree,
        }
        cycle.update(model_stats(refmac.xyzout))
        print("\nResidues built: %d" % cycle["residues_built"])
        print("Residues sequenced: %d" % cycle["residues_sequenced"])
        print("R-work: %.3f" % cycle["rwork"])
        print("R-free: %.3f" % cycle["rfree"])

        if self.improved(cycle):
            self.cycles_without_improvement = 0
        else:
            self.cycles_without_improvement += 1
            print("\nNo significant improvement for %d cycle(s)" % self.cycles_without_improvement)

        if cycle["rwork"] < self.min_rwork:
            self.min_rwork = cycle["rwork"]
        if cycle["rfree"] < self.min_rfree:
            self.min_rfree = cycle["rfree"]
            print("Copying files to output as R-free improved")
            shutil.copyfile(str(refmac.xyzout), "xyzout.pdb")
            shutil.copyfile(str(refmac.hklout.path), "hklout.mtz")
        if cycle["residues_built"] > self.max_residues_built:
            self.max_residues_built = cycle["residues_built"]
        if cycle["residues_sequenced"] > self.max_residues_sequenced:
            self.max_residues_sequenced = cycle["residues_sequenced"]
        if cycle["fragments_built"] < self.min_fragments_built:
            self.min_fragments_built = cycle["fragments_built"]
        if cycle["longest_fragment"] > self.max_longest_fragment:
            self.max_longest_fragment = cycle["longest_fragment"]

    def remove_job_directories(self, cycle):
        if self._args.keep_intermediate_files:
            return
        for job in self.jobs[cycle]:
            shutil.rmtree(job.directory)
