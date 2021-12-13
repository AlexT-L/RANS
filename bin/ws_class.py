from abc import ABC, abstractmethod

class WorkspaceClass(ABC):

    # check if dictionary has been initialized
    def __check_vars(self, workspace):
        if not workspace.has_dict(self.className):
            self.__init_vars(workspace)
