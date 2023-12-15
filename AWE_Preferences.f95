!
! AWE_Prefences.f95
!
! This file provides a set of functions to modify the default 
! behavior of an AWE application. These functions are called 
! before the primary AWE window opens.
!
! AWE_getStackSize
!
! This function returns the stacksize to use for an AWE application
INTEGER FUNCTION AWE_getStackSize() BIND(C, NAME="_AWE_getStackSize")
  IMPLICIT NONE
  
  AWE_getStackSize = (32 * 1024 * 1024)
  
  RETURN
END FUNCTION AWE_getStackSize

!
! AWE_getMdiMode
!
! This function controls whether windows opened in AWE will appear
! inside a single "frame" window or whether they open as individual
! windows.
!
LOGICAL FUNCTION AWE_getMdiMode() BIND(C, NAME="_AWE_getMdiMode")
  IMPLICIT NONE
  
  AWE_getMdiMode = .false.
  
  RETURN
END FUNCTION AWE_getMdiMode

!
! AWE_getShowMaximized
!
! This function can be used to open the AWE window already maximized.
!
LOGICAL FUNCTION AWE_getShowMaximized() BIND(C, NAME="_AWE_getShowMaximized")
  IMPLICIT NONE

  AWE_getShowMaximized = .false.

  RETURN
END FUNCTION AWE_getShowMaximized

!
! AWE_promptSaveOnExit
!
! This function controls whether AWE prompts to save the output
! window(s) at program exit. If this prompt is disabled, the contents
! of the window(s) will be lost if not explicitly saved.
!
LOGICAL FUNCTION AWE_promptSaveOnExit() BIND(C, NAME="_AWE_promptSaveOnExit")
  IMPLICIT NONE
  
  AWE_promptSaveOnExit = .false.
        
  RETURN
END FUNCTION AWE_promptSaveOnExit

!
! AWE_getMainWindowWidth
!
! This function controls the initial width of the window.
!
INTEGER FUNCTION AWE_getMainWindowWidth() BIND(C, NAME="_AWE_getMainWindowWidth")
  IMPLICIT NONE

  AWE_getMainWindowWidth = 1024
        
  RETURN
END FUNCTION AWE_getMainWindowWidth

!
! AWE_getMainWindowHeight
!
! This function controls the initial height of the window.
!
INTEGER FUNCTION AWE_getMainWindowHeight() BIND(C, NAME="_AWE_getMainWindowHeight")
  IMPLICIT NONE
  
   AWE_getMainWindowHeight = 768
        
   RETURN
END FUNCTION AWE_getMainWindowHeight

!
! AWE_defaultFontSize
!
! This function controls the height of the font use in the window.
!
INTEGER FUNCTION AWE_defaultFontSize() BIND(C, NAME="_AWE_defaultFontSize")
   IMPLICIT NONE

   AWE_defaultFontSize = 14

   RETURN
END FUNCTION AWE_defaultFontSize

!
! AWE_defaultFontFamily
!
! This function controls the family of the font use in the window.
!
SUBROUTINE AWE_defaultFontFamily(family) BIND(C, NAME="_AWE_defaultFontFamily")
   IMPLICIT NONE
   CHARACTER(LEN=*)  :: family
   
   family = "Courier New"
   
   RETURN
END SUBROUTINE AWE_defaultFontFamily

!
! _AWE_autoSave
!
! This function controls whether the window text is automatically
! saved when the program exits.
!
LOGICAL FUNCTION AWE_autoSave() BIND(C, NAME="_AWE_autoSave")
   IMPLICIT NONE
   
   AWE_autoSave = .false.
   
   RETURN
END FUNCTION AWE_autoSave

!
! _AWE_closeOnProgramFinish
!
! This function controls whether the AWE program exits immediately
! when the FORTAN program exits.  If false, the AWE windows will
! stay open until closed by the user.
!
LOGICAL FUNCTION AWE_closeOnProgramFinish() BIND(C, NAME="_AWE_closeOnProgramFinish")
   IMPLICIT NONE

   AWE_closeOnProgramFinish = .false.
        
   RETURN
END FUNCTION AWE_closeOnProgramFinish


!
! _AWE_showDefaultOutputWindow
!
! This function controls whether default AWE window is shown.  If
! you only want to show a plot or a canvas without the default text
! window, set this to .false.
!
LOGICAL FUNCTION AWE_showDefaultOutputWindow(), BIND(C, NAME="_AWE_showDefaultOutputWindow")
   IMPLICIT NONE

   AWE_showDefaultOutputWindow= .true.

   RETURN
END FUNCTION AWE_showDefaultOutputWindow

