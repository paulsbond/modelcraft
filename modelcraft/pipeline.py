from modelcraft.arguments import parse
from modelcraft.buccaneer import Buccaneer
from modelcraft.findwaters import FindWaters
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
        self.cycle = 0
        self.jobs = {0: []}
        self.current_hkl = self._args.hklin
        self.current_xyz = self._args.xyzin
        self.min_rwork = 1
        self.min_rfree = 1
        self.min_fragments_built = 999
        self.max_longest_fragment = 0
        self.max_residues_built = 0
        self.max_residues_sequenced = 0
        self.cycles_without_improvement = 0

    def _run(self):
        args = self._args
        if args.colin_hl is None and args.colin_phifom is None and args.mr_model is not None:
            print("\n## Preparations\n")
            self._get_phases_from_mr_model()
        for self.cycle in range(1, args.cycles + 1):
            print("\n## Cycle %d\n" % self.cycle)
            self.jobs[self.cycle] = []
            refmac = self._run_cycle()
            self.process_cycle_output(refmac)
            self.remove_job_directories(self.cycle - 1)
            if args.auto_stop and self.cycles_without_improvement == 4:
                break
        self.remove_job_directories(self.cycle)

    def _run_cycle(self):
        if self.cycle > 1:  # And resolution < 2.3 A?
            self._prune()
            self._refmac()
        self._buccaneer()
        self._refmac()
        self._prune(chains_only=True)
        if self._args.add_waters and self.min_rwork < 0.4:
            self._findwaters()
        return self._refmac()

    def _job_directory(self, name):
        directory = "%02d.%02d_%s" % (self.cycle, len(self.jobs[self.cycle]) + 1, name)
        print(directory)
        return directory

    def _get_phases_from_mr_model(self):
        directory = self._job_directory("mr_refinement")
        job = Refmac(self._args, directory, self._args.mr_model)
        self.jobs[self.cycle].append(job)
        self.current_hkl = job.hklout
        return job

    def _buccaneer(self):
        directory = self._job_directory("buccaneer")
        cycles = 3 if self.cycle == 1 else 2
        job = Buccaneer(self._args, directory, self.current_hkl, self.current_xyz, cycles)
        self.jobs[self.cycle].append(job)
        self.current_xyz = job.xyzout
        return job

    def _refmac(self):
        directory = self._job_directory("refmac")
        use_phases = self._args.unbiased and self.min_rwork > 0.35
        job = Refmac(self._args, directory, self.current_xyz, use_phases)
        self.jobs[self.cycle].append(job)
        self.current_hkl = job.hklout
        self.current_xyz = job.xyzout
        return job

    def _prune(self, chains_only=False):
        directory = self._job_directory("prune_chains" if chains_only else "prune")
        job = Prune(directory, self.current_xyz, self.current_hkl, chains_only)
        self.jobs[self.cycle].append(job)
        self.current_xyz = job.xyzout
        return job

    def _findwaters(self):
        directory = self._job_directory("findwaters")
        job = FindWaters(directory, self.current_xyz, self.current_hkl)
        self.jobs[self.cycle].append(job)
        self.current_xyz = job.xyzout
        return job

    def improved(self, refmac):
        required_improvement = 0.02
        improvement = (self.min_rwork - refmac.final_rwork) / self.min_rwork
        if improvement > required_improvement:
            return True
        improvement = (refmac.xyzout.residues - self.max_residues_built) / float(self.max_residues_built)
        if improvement > required_improvement:
            return True
        improvement = (refmac.xyzout.sequenced_residues - self.max_residues_sequenced) / float(self.max_residues_sequenced)
        if improvement > required_improvement:
            return True
        improvement = (self.min_fragments_built - refmac.xyzout.fragments) / float(self.min_fragments_built)
        if improvement > required_improvement:
            return True
        improvement = (refmac.xyzout.longest_fragment - self.max_longest_fragment) / float(self.max_longest_fragment)
        if improvement > required_improvement:
            return True
        return False

    def process_cycle_output(self, refmac):
        print("\nResidues built: %d" % refmac.xyzout.residues)
        print("Residues sequenced: %d" % refmac.xyzout.sequenced_residues)
        print("R-work: %.3f" % refmac.final_rwork)
        print("R-free: %.3f" % refmac.final_rfree)

        if self.improved(refmac):
            self.cycles_without_improvement = 0
        else:
            self.cycles_without_improvement += 1
            print("\nNo significant improvement for %d cycle(s)" % self.cycles_without_improvement)

        if refmac.final_rfree < self.min_rfree:
            self.min_rfree = refmac.final_rfree
            print("Copying files to output because R-free has improved")
            shutil.copyfile(str(refmac.xyzout.path), "xyzout.pdb")
            shutil.copyfile(str(refmac.hklout.path), "hklout.mtz")
        self.min_rwork = min(self.min_rwork, refmac.final_rwork)
        self.max_residues_built = max(self.max_residues_built, refmac.xyzout.residues)
        self.max_residues_sequenced = max(self.max_residues_sequenced, refmac.xyzout.sequenced_residues)
        self.min_fragments_built = min(self.min_fragments_built, refmac.xyzout.fragments)
        self.max_longest_fragment = max(self.max_longest_fragment, refmac.xyzout.longest_fragment)

    def remove_job_directories(self, cycle):
        if self._args.keep_intermediate_files:
            return
        for job in self.jobs[cycle]:
            shutil.rmtree(job.directory)
