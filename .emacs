;; Added by Package.el.  This must come before configurations of
;; installed packages.  Don't delete this line.  If you don't want it,
;; just comment it out by adding a semicolon to the start of the line.
;; You may delete these explanatory comments.
(package-initialize)			

(custom-set-variables
 ;; custom-set-variables was added by Custom.
 ;; If you edit it by hand, you could mess it up, so be careful.
 ;; Your init file should contain only one such instance.
 ;; If there is more than one, they won't work right.
 '(ansi-color-faces-vector
   [default default default italic underline success warning error])
 '(ansi-color-names-vector
   ["black" "red3" "ForestGreen" "yellow3" "blue" "magenta3" "DeepSkyBlue" "gray50"])
 '(custom-enabled-themes (quote (sanityinc-solarized-dark)))
 '(custom-safe-themes
   (quote
    ("4aee8551b53a43a883cb0b7f3255d6859d766b6c5e14bcb01bed572fcbef4328" default)))
 '(inhibit-startup-screen t)
 '(package-selected-packages (quote (smex ## company))))
(custom-set-faces
 ;; custom-set-faces was added by Custom.
 ;; If you edit it by hand, you could mess it up, so be careful.
 ;; Your init file should contain only one such instance.
 ;; If there is more than one, they won't work right.
 )

 ;;Line Number
(global-linum-mode t)
(column-number-mode t)
 
 ;;No Menu Bar
 ;(menu-bar-mode -1)
 ;;No Tool Bar
(tool-bar-mode -1)
 ;;No Scroll Bar
(scroll-bar-mode -1)
 ;;设置让光标指到某个括号与它匹配的括号
(show-paren-mode t)
(setq show-paren-style 'parentheses)

;(mouse-avoidance-mode 'animate)
;;光标靠近鼠标指针时，让鼠标指针自动让开，别挡住视线。

(setq frame-title-format "emacs@%b")
;;在标题栏显示buffer的名字，而不是 emacs@wangyin.com 这样没用的提示。
 ;;'y'/'n'代替'yes'/'no'
(fset 'yes-or-no-p 'y-or-n-p)

 ;;set transparent effect  
(global-set-key [(f8)] 'loop-alpha)  
(setq alpha-list '((100 100) (65 35)))  
(defun loop-alpha ()  
  (interactive)  
  (let ((h (car alpha-list)))                ;; head value will set to  
    ((lambda (a ab)  
       (set-frame-parameter (selected-frame) 'alpha (list a ab))  
       (add-to-list 'default-frame-alist (cons 'alpha (list a ab)))  
       ) (car h) (car (cdr h)))  
    (setq alpha-list (cdr (append alpha-list (list h))))  
  )  
)  
 ;; adding a package source
(setq package-archives '(("gnu" . "https://elpa.gnu.org/packages/")
                         ("marmalade" . "https://marmalade-repo.org/packages/")
                         ("melpa" . "https://melpa.org/packages/")))

 ;; company-mode
(add-hook 'after-init-hook 'global-company-mode) 

 ;; Evil
(add-to-list 'load-path "~/.emacs.d/evil")
(require 'evil)
(evil-mode 1)

 ;; Smex blind some keys
(global-set-key (kbd "M-x") 'smex)
(global-set-key (kbd "M-X") 'smex-major-mode-commands)
 ;; This is your old M-x.
(global-set-key (kbd "C-c C-c M-x") 'execute-extended-command)
