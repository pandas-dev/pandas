"use strict";
(self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] = self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] || []).push([[2816],{

/***/ 92816:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

// ESM COMPAT FLAG
__webpack_require__.r(__webpack_exports__);

// EXPORTS
__webpack_require__.d(__webpack_exports__, {
  Accordion: () => (/* reexport */ Accordion),
  AccordionItem: () => (/* reexport */ AccordionItem),
  Anchor: () => (/* reexport */ Anchor),
  AnchoredRegion: () => (/* reexport */ AnchoredRegion),
  Avatar: () => (/* reexport */ Avatar),
  Badge: () => (/* reexport */ Badge),
  Breadcrumb: () => (/* reexport */ Breadcrumb),
  BreadcrumbItem: () => (/* reexport */ BreadcrumbItem),
  Button: () => (/* reexport */ Button),
  Card: () => (/* reexport */ Card),
  Checkbox: () => (/* reexport */ Checkbox),
  Combobox: () => (/* reexport */ Combobox),
  DataGrid: () => (/* reexport */ DataGrid),
  DataGridCell: () => (/* reexport */ DataGridCell),
  DataGridRow: () => (/* reexport */ DataGridRow),
  DateField: () => (/* reexport */ DateField),
  Dialog: () => (/* reexport */ Dialog),
  Disclosure: () => (/* reexport */ Disclosure),
  Divider: () => (/* reexport */ Divider),
  Listbox: () => (/* reexport */ Listbox),
  Menu: () => (/* reexport */ Menu),
  MenuItem: () => (/* reexport */ MenuItem),
  NumberField: () => (/* reexport */ NumberField),
  Option: () => (/* reexport */ Option),
  Picker: () => (/* reexport */ Picker),
  PickerList: () => (/* reexport */ PickerList),
  PickerListItem: () => (/* reexport */ PickerListItem),
  PickerMenu: () => (/* reexport */ PickerMenu),
  PickerMenuOption: () => (/* reexport */ PickerMenuOption),
  Progress: () => (/* reexport */ Progress),
  ProgressRing: () => (/* reexport */ ProgressRing),
  Radio: () => (/* reexport */ Radio),
  RadioGroup: () => (/* reexport */ RadioGroup),
  Search: () => (/* reexport */ Search),
  Select: () => (/* reexport */ Select),
  Skeleton: () => (/* reexport */ Skeleton),
  Slider: () => (/* reexport */ Slider),
  SliderLabel: () => (/* reexport */ SliderLabel),
  Switch: () => (/* reexport */ Switch),
  Tab: () => (/* reexport */ Tab),
  TabPanel: () => (/* reexport */ TabPanel),
  Tabs: () => (/* reexport */ Tabs),
  TextArea: () => (/* reexport */ TextArea),
  TextField: () => (/* reexport */ TextField),
  Toolbar: () => (/* reexport */ Toolbar),
  Tooltip: () => (/* reexport */ Tooltip),
  TreeItem: () => (/* reexport */ TreeItem),
  TreeView: () => (/* reexport */ TreeView)
});

// EXTERNAL MODULE: consume shared module (default) @jupyter/web-components@~0.16.7 (singleton) (fallback: ../node_modules/@jupyter/web-components/dist/esm/index.js)
var index_js_ = __webpack_require__(83074);
// EXTERNAL MODULE: consume shared module (default) react@~18.2.0 (singleton) (fallback: ../node_modules/react/index.js)
var react_index_js_ = __webpack_require__(78156);
var react_index_js_default = /*#__PURE__*/__webpack_require__.n(react_index_js_);
;// CONCATENATED MODULE: ../node_modules/@jupyter/react-components/lib/react-utils.js


function useProperties(targetElement, propName, value) {
  (0,react_index_js_.useEffect)(() => {
    if (
      value !== undefined &&
      targetElement.current &&
      targetElement.current[propName] !== value
    ) {
      // add try catch to avoid errors when setting read-only properties
      try {
        targetElement.current[propName] = value;
      } catch (e) {
        console.warn(e);
      }
    }
  }, [value, targetElement.current]);
}

function useEventListener(targetElement, eventName, eventHandler) {
  (0,react_index_js_.useLayoutEffect)(() => {
    if (eventHandler !== undefined) {
      targetElement?.current?.addEventListener(eventName, eventHandler);
    }

    return () => {
      if (eventHandler?.cancel) {
        eventHandler.cancel();
      }

      targetElement?.current?.removeEventListener(eventName, eventHandler);
    };
  }, [eventName, eventHandler, targetElement.current]);
}

;// CONCATENATED MODULE: ../node_modules/@jupyter/react-components/lib/Accordion.js



(0,index_js_.provideJupyterDesignSystem)().register((0,index_js_.jpAccordion)());

const Accordion = (0,react_index_js_.forwardRef)((props, forwardedRef) => {
  const ref = (0,react_index_js_.useRef)(null);
  const { className, expandMode, ...filteredProps } = props;

  /** Event listeners - run once */
  useEventListener(ref, 'change', props.onChange);

  /** Properties - run whenever a property has changed */

  /** Methods - uses `useImperativeHandle` hook to pass ref to component */
  (0,react_index_js_.useImperativeHandle)(forwardedRef, () => ref.current, [ref.current]);

  return react_index_js_default().createElement(
    'jp-accordion',
    {
      ref,
      ...filteredProps,
      'expand-mode': props.expandMode || props['expand-mode'],
      class: props.className,
      exportparts: props.exportparts,
      for: props.htmlFor,
      part: props.part,
      tabindex: props.tabIndex,
      style: { ...props.style }
    },
    props.children
  );
});

;// CONCATENATED MODULE: ../node_modules/@jupyter/react-components/lib/AccordionItem.js



(0,index_js_.provideJupyterDesignSystem)().register((0,index_js_.jpAccordionItem)());

const AccordionItem = (0,react_index_js_.forwardRef)((props, forwardedRef) => {
  const ref = (0,react_index_js_.useRef)(null);
  const { className, headingLevel, id, expanded, ...filteredProps } = props;

  /** Event listeners - run once */
  useEventListener(ref, 'change', props.onChange);

  /** Properties - run whenever a property has changed */
  useProperties(ref, 'expanded', props.expanded);

  /** Methods - uses `useImperativeHandle` hook to pass ref to component */
  (0,react_index_js_.useImperativeHandle)(forwardedRef, () => ref.current, [ref.current]);

  return react_index_js_default().createElement(
    'jp-accordion-item',
    {
      ref,
      ...filteredProps,
      'heading-level': props.headingLevel || props['heading-level'],
      id: props.id,
      class: props.className,
      exportparts: props.exportparts,
      for: props.htmlFor,
      part: props.part,
      tabindex: props.tabIndex,
      style: { ...props.style }
    },
    props.children
  );
});

;// CONCATENATED MODULE: ../node_modules/@jupyter/react-components/lib/AnchoredRegion.js



(0,index_js_.provideJupyterDesignSystem)().register((0,index_js_.jpAnchoredRegion)());

const AnchoredRegion = (0,react_index_js_.forwardRef)((props, forwardedRef) => {
  const ref = (0,react_index_js_.useRef)(null);
  const {
    className,
    horizontalViewportLock,
    horizontalInset,
    verticalViewportLock,
    verticalInset,
    fixedPlacement,
    anchor,
    viewport,
    horizontalPositioningMode,
    horizontalDefaultPosition,
    horizontalThreshold,
    horizontalScaling,
    verticalPositioningMode,
    verticalDefaultPosition,
    verticalThreshold,
    verticalScaling,
    autoUpdateMode,
    anchorElement,
    viewportElement,
    verticalPosition,
    horizontalPosition,
    update,
    ...filteredProps
  } = props;

  /** Event listeners - run once */
  useEventListener(ref, 'loaded', props.onLoaded);
  useEventListener(ref, 'positionchange', props.onPositionchange);

  /** Properties - run whenever a property has changed */
  useProperties(ref, 'anchorElement', props.anchorElement);
  useProperties(ref, 'viewportElement', props.viewportElement);
  useProperties(ref, 'verticalPosition', props.verticalPosition);
  useProperties(ref, 'horizontalPosition', props.horizontalPosition);
  useProperties(ref, 'update', props.update);

  /** Methods - uses `useImperativeHandle` hook to pass ref to component */
  (0,react_index_js_.useImperativeHandle)(forwardedRef, () => ref.current, [ref.current]);

  return react_index_js_default().createElement(
    'jp-anchored-region',
    {
      ref,
      ...filteredProps,
      anchor: props.anchor,
      viewport: props.viewport,
      'horizontal-positioning-mode':
        props.horizontalPositioningMode || props['horizontal-positioning-mode'],
      'horizontal-default-position':
        props.horizontalDefaultPosition || props['horizontal-default-position'],
      'horizontal-threshold':
        props.horizontalThreshold || props['horizontal-threshold'],
      'horizontal-scaling':
        props.horizontalScaling || props['horizontal-scaling'],
      'vertical-positioning-mode':
        props.verticalPositioningMode || props['vertical-positioning-mode'],
      'vertical-default-position':
        props.verticalDefaultPosition || props['vertical-default-position'],
      'vertical-threshold':
        props.verticalThreshold || props['vertical-threshold'],
      'vertical-scaling': props.verticalScaling || props['vertical-scaling'],
      'auto-update-mode': props.autoUpdateMode || props['auto-update-mode'],
      class: props.className,
      exportparts: props.exportparts,
      for: props.htmlFor,
      part: props.part,
      tabindex: props.tabIndex,
      'horizontal-viewport-lock': props.horizontalViewportLock ? '' : undefined,
      'horizontal-inset': props.horizontalInset ? '' : undefined,
      'vertical-viewport-lock': props.verticalViewportLock ? '' : undefined,
      'vertical-inset': props.verticalInset ? '' : undefined,
      'fixed-placement': props.fixedPlacement ? '' : undefined,
      style: { ...props.style }
    },
    props.children
  );
});

;// CONCATENATED MODULE: ../node_modules/@jupyter/react-components/lib/Anchor.js



(0,index_js_.provideJupyterDesignSystem)().register((0,index_js_.jpAnchor)());

const Anchor = (0,react_index_js_.forwardRef)((props, forwardedRef) => {
  const ref = (0,react_index_js_.useRef)(null);
  const {
    className,
    appearance,
    download,
    href,
    hreflang,
    ping,
    referrerpolicy,
    rel,
    target,
    type,
    control,
    ...filteredProps
  } = props;

  /** Properties - run whenever a property has changed */
  useProperties(ref, 'control', props.control);

  /** Methods - uses `useImperativeHandle` hook to pass ref to component */
  (0,react_index_js_.useImperativeHandle)(forwardedRef, () => ref.current, [ref.current]);

  return react_index_js_default().createElement(
    'jp-anchor',
    {
      ref,
      ...filteredProps,
      appearance: props.appearance,
      download: props.download,
      href: props.href,
      hreflang: props.hreflang,
      ping: props.ping,
      referrerpolicy: props.referrerpolicy,
      rel: props.rel,
      target: props.target,
      type: props.type,
      class: props.className,
      exportparts: props.exportparts,
      for: props.htmlFor,
      part: props.part,
      tabindex: props.tabIndex,
      style: { ...props.style }
    },
    props.children
  );
});

;// CONCATENATED MODULE: ../node_modules/@jupyter/react-components/lib/Avatar.js



(0,index_js_.provideJupyterDesignSystem)().register((0,index_js_.jpAvatar)());

const Avatar = (0,react_index_js_.forwardRef)((props, forwardedRef) => {
  const ref = (0,react_index_js_.useRef)(null);
  const { className, src, alt, fill, color, link, shape, ...filteredProps } =
    props;

  /** Properties - run whenever a property has changed */

  /** Methods - uses `useImperativeHandle` hook to pass ref to component */
  (0,react_index_js_.useImperativeHandle)(forwardedRef, () => ref.current, [ref.current]);

  return react_index_js_default().createElement(
    'jp-avatar',
    {
      ref,
      ...filteredProps,
      src: props.src,
      alt: props.alt,
      fill: props.fill,
      color: props.color,
      link: props.link,
      shape: props.shape,
      class: props.className,
      exportparts: props.exportparts,
      for: props.htmlFor,
      part: props.part,
      tabindex: props.tabIndex,
      style: { ...props.style }
    },
    props.children
  );
});

;// CONCATENATED MODULE: ../node_modules/@jupyter/react-components/lib/Badge.js



(0,index_js_.provideJupyterDesignSystem)().register((0,index_js_.jpBadge)());

const Badge = (0,react_index_js_.forwardRef)((props, forwardedRef) => {
  const ref = (0,react_index_js_.useRef)(null);
  const { className, fill, color, circular, ...filteredProps } = props;

  /** Properties - run whenever a property has changed */
  useProperties(ref, 'circular', props.circular);

  /** Methods - uses `useImperativeHandle` hook to pass ref to component */
  (0,react_index_js_.useImperativeHandle)(forwardedRef, () => ref.current, [ref.current]);

  return react_index_js_default().createElement(
    'jp-badge',
    {
      ref,
      ...filteredProps,
      fill: props.fill,
      color: props.color,
      class: props.className,
      exportparts: props.exportparts,
      for: props.htmlFor,
      part: props.part,
      tabindex: props.tabIndex,
      style: { ...props.style }
    },
    props.children
  );
});

;// CONCATENATED MODULE: ../node_modules/@jupyter/react-components/lib/Breadcrumb.js



(0,index_js_.provideJupyterDesignSystem)().register((0,index_js_.jpBreadcrumb)());

const Breadcrumb = (0,react_index_js_.forwardRef)((props, forwardedRef) => {
  const ref = (0,react_index_js_.useRef)(null);
  const { className, ...filteredProps } = props;

  /** Properties - run whenever a property has changed */

  /** Methods - uses `useImperativeHandle` hook to pass ref to component */
  (0,react_index_js_.useImperativeHandle)(forwardedRef, () => ref.current, [ref.current]);

  return react_index_js_default().createElement(
    'jp-breadcrumb',
    {
      ref,
      ...filteredProps,
      class: props.className,
      exportparts: props.exportparts,
      for: props.htmlFor,
      part: props.part,
      tabindex: props.tabIndex,
      style: { ...props.style }
    },
    props.children
  );
});

;// CONCATENATED MODULE: ../node_modules/@jupyter/react-components/lib/BreadcrumbItem.js



(0,index_js_.provideJupyterDesignSystem)().register((0,index_js_.jpBreadcrumbItem)());

const BreadcrumbItem = (0,react_index_js_.forwardRef)((props, forwardedRef) => {
  const ref = (0,react_index_js_.useRef)(null);
  const {
    className,
    download,
    href,
    hreflang,
    ping,
    referrerpolicy,
    rel,
    target,
    type,
    control,
    ...filteredProps
  } = props;

  /** Properties - run whenever a property has changed */
  useProperties(ref, 'control', props.control);

  /** Methods - uses `useImperativeHandle` hook to pass ref to component */
  (0,react_index_js_.useImperativeHandle)(forwardedRef, () => ref.current, [ref.current]);

  return react_index_js_default().createElement(
    'jp-breadcrumb-item',
    {
      ref,
      ...filteredProps,
      download: props.download,
      href: props.href,
      hreflang: props.hreflang,
      ping: props.ping,
      referrerpolicy: props.referrerpolicy,
      rel: props.rel,
      target: props.target,
      type: props.type,
      class: props.className,
      exportparts: props.exportparts,
      for: props.htmlFor,
      part: props.part,
      tabindex: props.tabIndex,
      style: { ...props.style }
    },
    props.children
  );
});

;// CONCATENATED MODULE: ../node_modules/@jupyter/react-components/lib/Button.js



(0,index_js_.provideJupyterDesignSystem)().register((0,index_js_.jpButton)());

const Button = (0,react_index_js_.forwardRef)((props, forwardedRef) => {
  const ref = (0,react_index_js_.useRef)(null);
  const {
    className,
    minimal,
    appearance,
    form,
    formaction,
    formenctype,
    formmethod,
    formtarget,
    type,
    autofocus,
    formnovalidate,
    defaultSlottedContent,
    disabled,
    required,
    ...filteredProps
  } = props;

  /** Properties - run whenever a property has changed */
  useProperties(ref, 'autofocus', props.autofocus);
  useProperties(ref, 'formnovalidate', props.formnovalidate);
  useProperties(ref, 'defaultSlottedContent', props.defaultSlottedContent);
  useProperties(ref, 'disabled', props.disabled);
  useProperties(ref, 'required', props.required);

  /** Methods - uses `useImperativeHandle` hook to pass ref to component */
  (0,react_index_js_.useImperativeHandle)(forwardedRef, () => ref.current, [ref.current]);

  return react_index_js_default().createElement(
    'jp-button',
    {
      ref,
      ...filteredProps,
      appearance: props.appearance,
      form: props.form,
      formaction: props.formaction,
      formenctype: props.formenctype,
      formmethod: props.formmethod,
      formtarget: props.formtarget,
      type: props.type,
      class: props.className,
      exportparts: props.exportparts,
      for: props.htmlFor,
      part: props.part,
      tabindex: props.tabIndex,
      minimal: props.minimal ? '' : undefined,
      style: { ...props.style }
    },
    props.children
  );
});

;// CONCATENATED MODULE: ../node_modules/@jupyter/react-components/lib/Card.js



(0,index_js_.provideJupyterDesignSystem)().register((0,index_js_.jpCard)());

const Card = (0,react_index_js_.forwardRef)((props, forwardedRef) => {
  const ref = (0,react_index_js_.useRef)(null);
  const { className, ...filteredProps } = props;

  /** Properties - run whenever a property has changed */

  /** Methods - uses `useImperativeHandle` hook to pass ref to component */
  (0,react_index_js_.useImperativeHandle)(forwardedRef, () => ref.current, [ref.current]);

  return react_index_js_default().createElement(
    'jp-card',
    {
      ref,
      ...filteredProps,
      class: props.className,
      exportparts: props.exportparts,
      for: props.htmlFor,
      part: props.part,
      tabindex: props.tabIndex,
      style: { ...props.style }
    },
    props.children
  );
});

;// CONCATENATED MODULE: ../node_modules/@jupyter/react-components/lib/Checkbox.js



(0,index_js_.provideJupyterDesignSystem)().register((0,index_js_.jpCheckbox)());

const Checkbox = (0,react_index_js_.forwardRef)((props, forwardedRef) => {
  const ref = (0,react_index_js_.useRef)(null);
  const {
    className,
    readonly,
    readOnly,
    indeterminate,
    checked,
    disabled,
    required,
    ...filteredProps
  } = props;

  /** Event listeners - run once */
  useEventListener(ref, 'change', props.onChange);

  /** Properties - run whenever a property has changed */
  useProperties(ref, 'readOnly', props.readOnly);
  useProperties(ref, 'indeterminate', props.indeterminate);
  useProperties(ref, 'checked', props.checked);
  useProperties(ref, 'disabled', props.disabled);
  useProperties(ref, 'required', props.required);

  /** Methods - uses `useImperativeHandle` hook to pass ref to component */
  (0,react_index_js_.useImperativeHandle)(forwardedRef, () => ref.current, [ref.current]);

  // Add web component internal classes on top of `className`
  let allClasses = className ?? '';
  if (ref.current?.indeterminate) {
    allClasses += ' indeterminate';
  }

  return react_index_js_default().createElement(
    'jp-checkbox',
    {
      ref,
      ...filteredProps,
      class: allClasses.trim(),
      exportparts: props.exportparts,
      for: props.htmlFor,
      part: props.part,
      tabindex: props.tabIndex,
      readonly: props.readonly ? '' : undefined,
      style: { ...props.style }
    },
    props.children
  );
});

;// CONCATENATED MODULE: ../node_modules/@jupyter/react-components/lib/Combobox.js



(0,index_js_.provideJupyterDesignSystem)().register((0,index_js_.jpCombobox)());

const Combobox = (0,react_index_js_.forwardRef)((props, forwardedRef) => {
  const ref = (0,react_index_js_.useRef)(null);
  const {
    className,
    autowidth,
    minimal,
    open,
    autocomplete,
    placeholder,
    position,
    autoWidth,
    filteredOptions,
    options,
    value,
    length,
    disabled,
    selectedIndex,
    selectedOptions,
    required,
    ...filteredProps
  } = props;

  /** Event listeners - run once */
  useEventListener(ref, 'input', props.onInput);
  useEventListener(ref, 'change', props.onChange);

  /** Properties - run whenever a property has changed */
  useProperties(ref, 'autoWidth', props.autoWidth);
  useProperties(ref, 'filteredOptions', props.filteredOptions);
  useProperties(ref, 'options', props.options);
  useProperties(ref, 'value', props.value);
  useProperties(ref, 'length', props.length);
  useProperties(ref, 'disabled', props.disabled);
  useProperties(ref, 'selectedIndex', props.selectedIndex);
  useProperties(ref, 'selectedOptions', props.selectedOptions);
  useProperties(ref, 'required', props.required);

  /** Methods - uses `useImperativeHandle` hook to pass ref to component */
  (0,react_index_js_.useImperativeHandle)(forwardedRef, () => ref.current, [ref.current]);

  return react_index_js_default().createElement(
    'jp-combobox',
    {
      ref,
      ...filteredProps,
      autocomplete: props.autocomplete,
      placeholder: props.placeholder,
      position: props.position,
      class: props.className,
      exportparts: props.exportparts,
      for: props.htmlFor,
      part: props.part,
      tabindex: props.tabIndex,
      autowidth: props.autowidth ? '' : undefined,
      minimal: props.minimal ? '' : undefined,
      open: props.open ? '' : undefined,
      style: { ...props.style }
    },
    props.children
  );
});

;// CONCATENATED MODULE: ../node_modules/@jupyter/react-components/lib/DateField.js



(0,index_js_.provideJupyterDesignSystem)().register((0,index_js_.jpDateField)());

const DateField = (0,react_index_js_.forwardRef)((props, forwardedRef) => {
  const ref = (0,react_index_js_.useRef)(null);
  const {
    className,
    autofocus,
    step,
    max,
    min,
    disabled,
    required,
    ...filteredProps
  } = props;

  /** Event listeners - run once */
  useEventListener(ref, 'input', props.onInput);
  useEventListener(ref, 'change', props.onChange);

  /** Properties - run whenever a property has changed */
  useProperties(ref, 'autofocus', props.autofocus);
  useProperties(ref, 'step', props.step);
  useProperties(ref, 'max', props.max);
  useProperties(ref, 'min', props.min);
  useProperties(ref, 'disabled', props.disabled);
  useProperties(ref, 'required', props.required);

  /** Methods - uses `useImperativeHandle` hook to pass ref to component */
  (0,react_index_js_.useImperativeHandle)(forwardedRef, () => ref.current, [ref.current]);

  return react_index_js_default().createElement(
    'jp-date-field',
    {
      ref,
      ...filteredProps,
      class: props.className,
      exportparts: props.exportparts,
      for: props.htmlFor,
      part: props.part,
      tabindex: props.tabIndex,
      style: { ...props.style }
    },
    props.children
  );
});

;// CONCATENATED MODULE: ../node_modules/@jupyter/react-components/lib/DataGridCell.js



(0,index_js_.provideJupyterDesignSystem)().register((0,index_js_.jpDataGridCell)());

const DataGridCell = (0,react_index_js_.forwardRef)((props, forwardedRef) => {
  const ref = (0,react_index_js_.useRef)(null);
  const {
    className,
    cellType,
    gridColumn,
    rowData,
    columnDefinition,
    ...filteredProps
  } = props;

  /** Event listeners - run once */
  useEventListener(ref, 'cell-focused', props.onCellFocused);

  /** Properties - run whenever a property has changed */
  useProperties(ref, 'rowData', props.rowData);
  useProperties(ref, 'columnDefinition', props.columnDefinition);

  /** Methods - uses `useImperativeHandle` hook to pass ref to component */
  (0,react_index_js_.useImperativeHandle)(forwardedRef, () => ref.current, [ref.current]);

  // Add web component internal classes on top of `className`
  let allClasses = className ?? '';
  if (ref.current?.cellType === 'columnheader') {
    allClasses += ' column-header';
  }

  return react_index_js_default().createElement(
    'jp-data-grid-cell',
    {
      ref,
      ...filteredProps,
      'cell-type': props.cellType || props['cell-type'],
      'grid-column': props.gridColumn || props['grid-column'],
      class: allClasses.trim(),
      exportparts: props.exportparts,
      for: props.htmlFor,
      part: props.part,
      tabindex: props.tabIndex,
      style: { ...props.style }
    },
    props.children
  );
});

;// CONCATENATED MODULE: ../node_modules/@jupyter/react-components/lib/DataGridRow.js



(0,index_js_.provideJupyterDesignSystem)().register((0,index_js_.jpDataGridRow)());

const DataGridRow = (0,react_index_js_.forwardRef)((props, forwardedRef) => {
  const ref = (0,react_index_js_.useRef)(null);
  const {
    className,
    gridTemplateColumns,
    rowType,
    rowData,
    columnDefinitions,
    cellItemTemplate,
    headerCellItemTemplate,
    rowIndex,
    ...filteredProps
  } = props;

  /** Event listeners - run once */
  useEventListener(ref, 'row-focused', props.onRowFocused);

  /** Properties - run whenever a property has changed */
  useProperties(ref, 'rowData', props.rowData);
  useProperties(ref, 'columnDefinitions', props.columnDefinitions);
  useProperties(ref, 'cellItemTemplate', props.cellItemTemplate);
  useProperties(ref, 'headerCellItemTemplate', props.headerCellItemTemplate);
  useProperties(ref, 'rowIndex', props.rowIndex);

  /** Methods - uses `useImperativeHandle` hook to pass ref to component */
  (0,react_index_js_.useImperativeHandle)(forwardedRef, () => ref.current, [ref.current]);

  // Add web component internal classes on top of `className`
  let allClasses = className ?? '';
  if (ref.current) {
    if (ref.current.rowType !== 'default') {
      allClasses += ` ${ref.current.rowType}`;
    }
  }

  return react_index_js_default().createElement(
    'jp-data-grid-row',
    {
      ref,
      ...filteredProps,
      'grid-template-columns':
        props.gridTemplateColumns || props['grid-template-columns'],
      'row-type': props.rowType || props['row-type'],
      class: allClasses.trim(),
      exportparts: props.exportparts,
      for: props.htmlFor,
      part: props.part,
      tabindex: props.tabIndex,
      style: { ...props.style }
    },
    props.children
  );
});

;// CONCATENATED MODULE: ../node_modules/@jupyter/react-components/lib/DataGrid.js



(0,index_js_.provideJupyterDesignSystem)().register((0,index_js_.jpDataGrid)());

const DataGrid = (0,react_index_js_.forwardRef)((props, forwardedRef) => {
  const ref = (0,react_index_js_.useRef)(null);
  const {
    className,
    noTabbing,
    generateHeader,
    gridTemplateColumns,
    rowsData,
    columnDefinitions,
    rowItemTemplate,
    cellItemTemplate,
    headerCellItemTemplate,
    focusRowIndex,
    focusColumnIndex,
    rowElementTag,
    ...filteredProps
  } = props;

  /** Properties - run whenever a property has changed */
  useProperties(ref, 'rowsData', props.rowsData);
  useProperties(ref, 'columnDefinitions', props.columnDefinitions);
  useProperties(ref, 'rowItemTemplate', props.rowItemTemplate);
  useProperties(ref, 'cellItemTemplate', props.cellItemTemplate);
  useProperties(ref, 'headerCellItemTemplate', props.headerCellItemTemplate);
  useProperties(ref, 'focusRowIndex', props.focusRowIndex);
  useProperties(ref, 'focusColumnIndex', props.focusColumnIndex);
  useProperties(ref, 'rowElementTag', props.rowElementTag);

  /** Methods - uses `useImperativeHandle` hook to pass ref to component */
  (0,react_index_js_.useImperativeHandle)(forwardedRef, () => ref.current, [ref.current]);

  return react_index_js_default().createElement(
    'jp-data-grid',
    {
      ref,
      ...filteredProps,
      'generate-header': props.generateHeader || props['generate-header'],
      'grid-template-columns':
        props.gridTemplateColumns || props['grid-template-columns'],
      class: props.className,
      exportparts: props.exportparts,
      for: props.htmlFor,
      part: props.part,
      tabindex: props.tabIndex,
      'no-tabbing': props.noTabbing ? '' : undefined,
      style: { ...props.style }
    },
    props.children
  );
});

;// CONCATENATED MODULE: ../node_modules/@jupyter/react-components/lib/Dialog.js



(0,index_js_.provideJupyterDesignSystem)().register((0,index_js_.jpDialog)());

const Dialog = (0,react_index_js_.forwardRef)((props, forwardedRef) => {
  const ref = (0,react_index_js_.useRef)(null);
  const {
    className,
    trapFocus,
    ariaDescribedby,
    ariaLabelledby,
    ariaLabel,
    modal,
    hidden,
    ...filteredProps
  } = props;

  /** Event listeners - run once */
  useEventListener(ref, 'cancel', props.onCancel);
  useEventListener(ref, 'close', props.onClose);

  /** Properties - run whenever a property has changed */
  useProperties(ref, 'modal', props.modal);
  useProperties(ref, 'hidden', props.hidden);

  /** Methods - uses `useImperativeHandle` hook to pass ref to component */
  (0,react_index_js_.useImperativeHandle)(forwardedRef, () => ({
    show: () => ref.current.show(),
    hide: () => ref.current.hide(),
    compose: (this_, elementDefinition) =>
      ref.current.compose(this_, elementDefinition)
  }));

  return react_index_js_default().createElement(
    'jp-dialog',
    {
      ref,
      ...filteredProps,
      'aria-describedby': props.ariaDescribedby || props['aria-describedby'],
      'aria-labelledby': props.ariaLabelledby || props['aria-labelledby'],
      'aria-label': props.ariaLabel || props['aria-label'],
      class: props.className,
      exportparts: props.exportparts,
      for: props.htmlFor,
      part: props.part,
      tabindex: props.tabIndex,
      'trap-focus': props.trapFocus ? '' : undefined,
      style: { ...props.style }
    },
    props.children
  );
});

;// CONCATENATED MODULE: ../node_modules/@jupyter/react-components/lib/Disclosure.js



(0,index_js_.provideJupyterDesignSystem)().register((0,index_js_.jpDisclosure)());

const Disclosure = (0,react_index_js_.forwardRef)((props, forwardedRef) => {
  const ref = (0,react_index_js_.useRef)(null);
  const { className, appearance, title, expanded, ...filteredProps } = props;

  /** Event listeners - run once */
  useEventListener(ref, 'toggle', props.onToggle);

  /** Properties - run whenever a property has changed */
  useProperties(ref, 'expanded', props.expanded);

  /** Methods - uses `useImperativeHandle` hook to pass ref to component */
  (0,react_index_js_.useImperativeHandle)(forwardedRef, () => ref.current, [ref.current]);

  return react_index_js_default().createElement(
    'jp-disclosure',
    {
      ref,
      ...filteredProps,
      appearance: props.appearance,
      title: props.title,
      class: props.className,
      exportparts: props.exportparts,
      for: props.htmlFor,
      part: props.part,
      tabindex: props.tabIndex,
      style: { ...props.style }
    },
    props.children
  );
});

;// CONCATENATED MODULE: ../node_modules/@jupyter/react-components/lib/Divider.js



(0,index_js_.provideJupyterDesignSystem)().register((0,index_js_.jpDivider)());

const Divider = (0,react_index_js_.forwardRef)((props, forwardedRef) => {
  const ref = (0,react_index_js_.useRef)(null);
  const { className, role, orientation, ...filteredProps } = props;

  /** Properties - run whenever a property has changed */

  /** Methods - uses `useImperativeHandle` hook to pass ref to component */
  (0,react_index_js_.useImperativeHandle)(forwardedRef, () => ref.current, [ref.current]);

  return react_index_js_default().createElement(
    'jp-divider',
    {
      ref,
      ...filteredProps,
      role: props.role,
      orientation: props.orientation,
      class: props.className,
      exportparts: props.exportparts,
      for: props.htmlFor,
      part: props.part,
      tabindex: props.tabIndex,
      style: { ...props.style }
    },
    props.children
  );
});

;// CONCATENATED MODULE: ../node_modules/@jupyter/react-components/lib/Listbox.js



(0,index_js_.provideJupyterDesignSystem)().register((0,index_js_.jpListbox)());

const Listbox = (0,react_index_js_.forwardRef)((props, forwardedRef) => {
  const ref = (0,react_index_js_.useRef)(null);
  const {
    className,
    multiple,
    size,
    length,
    options,
    disabled,
    selectedIndex,
    selectedOptions,
    ...filteredProps
  } = props;

  /** Event listeners - run once */
  useEventListener(ref, 'change', props.onChange);

  /** Properties - run whenever a property has changed */
  useProperties(ref, 'multiple', props.multiple);
  useProperties(ref, 'size', props.size);
  useProperties(ref, 'length', props.length);
  useProperties(ref, 'options', props.options);
  useProperties(ref, 'disabled', props.disabled);
  useProperties(ref, 'selectedIndex', props.selectedIndex);
  useProperties(ref, 'selectedOptions', props.selectedOptions);

  /** Methods - uses `useImperativeHandle` hook to pass ref to component */
  (0,react_index_js_.useImperativeHandle)(forwardedRef, () => ref.current, [ref.current]);

  return react_index_js_default().createElement(
    'jp-listbox',
    {
      ref,
      ...filteredProps,
      class: props.className,
      exportparts: props.exportparts,
      for: props.htmlFor,
      part: props.part,
      tabindex: props.tabIndex,
      style: { ...props.style }
    },
    props.children
  );
});

;// CONCATENATED MODULE: ../node_modules/@jupyter/react-components/lib/MenuItem.js



(0,index_js_.provideJupyterDesignSystem)().register((0,index_js_.jpMenuItem)());

const MenuItem = (0,react_index_js_.forwardRef)((props, forwardedRef) => {
  const ref = (0,react_index_js_.useRef)(null);
  const { className, role, disabled, expanded, checked, ...filteredProps } =
    props;

  /** Event listeners - run once */
  useEventListener(ref, 'expanded-change', props.onExpand);
  useEventListener(ref, 'change', props.onChange);

  /** Properties - run whenever a property has changed */
  useProperties(ref, 'disabled', props.disabled);
  useProperties(ref, 'expanded', props.expanded);
  useProperties(ref, 'checked', props.checked);

  /** Methods - uses `useImperativeHandle` hook to pass ref to component */
  (0,react_index_js_.useImperativeHandle)(forwardedRef, () => ref.current, [ref.current]);

  // Add web component internal classes on top of `className`
  let allClasses = className ?? '';
  if (ref.current) {
    allClasses += ` indent-${ref.current.startColumnCount}`;
    if (ref.current.expanded) {
      allClasses += ' expanded';
    }
  }

  return react_index_js_default().createElement(
    'jp-menu-item',
    {
      ref,
      ...filteredProps,
      role: props.role,
      class: allClasses.trim(),
      exportparts: props.exportparts,
      for: props.htmlFor,
      part: props.part,
      tabindex: props.tabIndex,
      style: { ...props.style }
    },
    props.children
  );
});

;// CONCATENATED MODULE: ../node_modules/@jupyter/react-components/lib/Menu.js



(0,index_js_.provideJupyterDesignSystem)().register((0,index_js_.jpMenu)());

const Menu = (0,react_index_js_.forwardRef)((props, forwardedRef) => {
  const ref = (0,react_index_js_.useRef)(null);
  const { className, ...filteredProps } = props;

  /** Properties - run whenever a property has changed */

  /** Methods - uses `useImperativeHandle` hook to pass ref to component */
  (0,react_index_js_.useImperativeHandle)(forwardedRef, () => ref.current, [ref.current]);

  return react_index_js_default().createElement(
    'jp-menu',
    {
      ref,
      ...filteredProps,
      class: props.className,
      exportparts: props.exportparts,
      for: props.htmlFor,
      part: props.part,
      tabindex: props.tabIndex,
      style: { ...props.style }
    },
    props.children
  );
});

;// CONCATENATED MODULE: ../node_modules/@jupyter/react-components/lib/NumberField.js



(0,index_js_.provideJupyterDesignSystem)().register((0,index_js_.jpNumberField)());

const NumberField = (0,react_index_js_.forwardRef)((props, forwardedRef) => {
  const ref = (0,react_index_js_.useRef)(null);
  const {
    className,
    readonly,
    hideStep,
    appearance,
    placeholder,
    list,
    readOnly,
    autofocus,
    maxlength,
    minlength,
    size,
    step,
    max,
    min,
    disabled,
    required,
    ...filteredProps
  } = props;

  /** Event listeners - run once */
  useEventListener(ref, 'input', props.onInput);
  useEventListener(ref, 'change', props.onChange);

  /** Properties - run whenever a property has changed */
  useProperties(ref, 'readOnly', props.readOnly);
  useProperties(ref, 'autofocus', props.autofocus);
  useProperties(ref, 'maxlength', props.maxlength);
  useProperties(ref, 'minlength', props.minlength);
  useProperties(ref, 'size', props.size);
  useProperties(ref, 'step', props.step);
  useProperties(ref, 'max', props.max);
  useProperties(ref, 'min', props.min);
  useProperties(ref, 'disabled', props.disabled);
  useProperties(ref, 'required', props.required);

  /** Methods - uses `useImperativeHandle` hook to pass ref to component */
  (0,react_index_js_.useImperativeHandle)(forwardedRef, () => ref.current, [ref.current]);

  return react_index_js_default().createElement(
    'jp-number-field',
    {
      ref,
      ...filteredProps,
      appearance: props.appearance,
      placeholder: props.placeholder,
      list: props.list,
      class: props.className,
      exportparts: props.exportparts,
      for: props.htmlFor,
      part: props.part,
      tabindex: props.tabIndex,
      readonly: props.readonly ? '' : undefined,
      'hide-step': props.hideStep ? '' : undefined,
      style: { ...props.style }
    },
    props.children
  );
});

;// CONCATENATED MODULE: ../node_modules/@jupyter/react-components/lib/Option.js



(0,index_js_.provideJupyterDesignSystem)().register((0,index_js_.jpOption)());

const Option = (0,react_index_js_.forwardRef)((props, forwardedRef) => {
  const ref = (0,react_index_js_.useRef)(null);
  const {
    className,
    selected,
    value,
    checked,
    content,
    defaultSelected,
    disabled,
    selectedAttribute,
    dirtyValue,
    ...filteredProps
  } = props;

  /** Properties - run whenever a property has changed */
  useProperties(ref, 'checked', props.checked);
  useProperties(ref, 'content', props.content);
  useProperties(ref, 'defaultSelected', props.defaultSelected);
  useProperties(ref, 'disabled', props.disabled);
  useProperties(ref, 'selectedAttribute', props.selectedAttribute);
  useProperties(ref, 'dirtyValue', props.dirtyValue);

  /** Methods - uses `useImperativeHandle` hook to pass ref to component */
  (0,react_index_js_.useImperativeHandle)(forwardedRef, () => ref.current, [ref.current]);

  return react_index_js_default().createElement(
    'jp-option',
    {
      ref,
      ...filteredProps,
      value: props.value,
      class: props.className,
      exportparts: props.exportparts,
      for: props.htmlFor,
      part: props.part,
      tabindex: props.tabIndex,
      selected: props.selected ? '' : undefined,
      style: { ...props.style }
    },
    props.children
  );
});

;// CONCATENATED MODULE: ../node_modules/@jupyter/react-components/lib/ProgressRing.js



(0,index_js_.provideJupyterDesignSystem)().register((0,index_js_.jpProgressRing)());

const ProgressRing = (0,react_index_js_.forwardRef)((props, forwardedRef) => {
  const ref = (0,react_index_js_.useRef)(null);
  const { className, value, min, max, paused, ...filteredProps } = props;

  /** Properties - run whenever a property has changed */
  useProperties(ref, 'value', props.value);
  useProperties(ref, 'min', props.min);
  useProperties(ref, 'max', props.max);
  useProperties(ref, 'paused', props.paused);

  /** Methods - uses `useImperativeHandle` hook to pass ref to component */
  (0,react_index_js_.useImperativeHandle)(forwardedRef, () => ref.current, [ref.current]);

  return react_index_js_default().createElement(
    'jp-progress-ring',
    {
      ref,
      ...filteredProps,
      class: props.className,
      exportparts: props.exportparts,
      for: props.htmlFor,
      part: props.part,
      tabindex: props.tabIndex,
      style: { ...props.style }
    },
    props.children
  );
});

;// CONCATENATED MODULE: ../node_modules/@jupyter/react-components/lib/Progress.js



(0,index_js_.provideJupyterDesignSystem)().register((0,index_js_.jpProgress)());

const Progress = (0,react_index_js_.forwardRef)((props, forwardedRef) => {
  const ref = (0,react_index_js_.useRef)(null);
  const { className, value, min, max, paused, ...filteredProps } = props;

  /** Properties - run whenever a property has changed */
  useProperties(ref, 'value', props.value);
  useProperties(ref, 'min', props.min);
  useProperties(ref, 'max', props.max);
  useProperties(ref, 'paused', props.paused);

  /** Methods - uses `useImperativeHandle` hook to pass ref to component */
  (0,react_index_js_.useImperativeHandle)(forwardedRef, () => ref.current, [ref.current]);

  return react_index_js_default().createElement(
    'jp-progress',
    {
      ref,
      ...filteredProps,
      class: props.className,
      exportparts: props.exportparts,
      for: props.htmlFor,
      part: props.part,
      tabindex: props.tabIndex,
      style: { ...props.style }
    },
    props.children
  );
});

;// CONCATENATED MODULE: ../node_modules/@jupyter/react-components/lib/Radio.js



(0,index_js_.provideJupyterDesignSystem)().register((0,index_js_.jpRadio)());

const Radio = (0,react_index_js_.forwardRef)((props, forwardedRef) => {
  const ref = (0,react_index_js_.useRef)(null);
  const {
    className,
    readonly,
    readOnly,
    name,
    checked,
    disabled,
    required,
    ...filteredProps
  } = props;

  /** Event listeners - run once */
  useEventListener(ref, 'change', props.onChange);

  /** Properties - run whenever a property has changed */
  useProperties(ref, 'readOnly', props.readOnly);
  useProperties(ref, 'name', props.name);
  useProperties(ref, 'checked', props.checked);
  useProperties(ref, 'disabled', props.disabled);
  useProperties(ref, 'required', props.required);

  /** Methods - uses `useImperativeHandle` hook to pass ref to component */
  (0,react_index_js_.useImperativeHandle)(forwardedRef, () => ref.current, [ref.current]);

  return react_index_js_default().createElement(
    'jp-radio',
    {
      ref,
      ...filteredProps,
      class: props.className,
      exportparts: props.exportparts,
      for: props.htmlFor,
      part: props.part,
      tabindex: props.tabIndex,
      readonly: props.readonly ? '' : undefined,
      style: { ...props.style }
    },
    props.children
  );
});

;// CONCATENATED MODULE: ../node_modules/@jupyter/react-components/lib/RadioGroup.js



(0,index_js_.provideJupyterDesignSystem)().register((0,index_js_.jpRadioGroup)());

const RadioGroup = (0,react_index_js_.forwardRef)((props, forwardedRef) => {
  const ref = (0,react_index_js_.useRef)(null);
  const {
    className,
    readonly,
    disabled,
    name,
    value,
    orientation,
    readOnly,
    ...filteredProps
  } = props;

  /** Event listeners - run once */
  useEventListener(ref, 'change', props.onChange);

  /** Properties - run whenever a property has changed */
  useProperties(ref, 'readOnly', props.readOnly);

  /** Methods - uses `useImperativeHandle` hook to pass ref to component */
  (0,react_index_js_.useImperativeHandle)(forwardedRef, () => ref.current, [ref.current]);

  return react_index_js_default().createElement(
    'jp-radio-group',
    {
      ref,
      ...filteredProps,
      name: props.name,
      value: props.value,
      orientation: props.orientation,
      class: props.className,
      exportparts: props.exportparts,
      for: props.htmlFor,
      part: props.part,
      tabindex: props.tabIndex,
      readonly: props.readonly ? '' : undefined,
      disabled: props.disabled ? '' : undefined,
      style: { ...props.style }
    },
    props.children
  );
});

;// CONCATENATED MODULE: ../node_modules/@jupyter/react-components/lib/Search.js



(0,index_js_.provideJupyterDesignSystem)().register((0,index_js_.jpSearch)());

const Search = (0,react_index_js_.forwardRef)((props, forwardedRef) => {
  const ref = (0,react_index_js_.useRef)(null);
  const {
    className,
    readonly,
    appearance,
    placeholder,
    list,
    pattern,
    readOnly,
    autofocus,
    maxlength,
    minlength,
    size,
    spellcheck,
    disabled,
    required,
    ...filteredProps
  } = props;

  /** Event listeners - run once */
  useEventListener(ref, 'input', props.onInput);
  useEventListener(ref, 'change', props.onChange);

  /** Properties - run whenever a property has changed */
  useProperties(ref, 'readOnly', props.readOnly);
  useProperties(ref, 'autofocus', props.autofocus);
  useProperties(ref, 'maxlength', props.maxlength);
  useProperties(ref, 'minlength', props.minlength);
  useProperties(ref, 'size', props.size);
  useProperties(ref, 'spellcheck', props.spellcheck);
  useProperties(ref, 'disabled', props.disabled);
  useProperties(ref, 'required', props.required);

  /** Methods - uses `useImperativeHandle` hook to pass ref to component */
  (0,react_index_js_.useImperativeHandle)(forwardedRef, () => ref.current, [ref.current]);

  return react_index_js_default().createElement(
    'jp-search',
    {
      ref,
      ...filteredProps,
      appearance: props.appearance,
      placeholder: props.placeholder,
      list: props.list,
      pattern: props.pattern,
      class: props.className,
      exportparts: props.exportparts,
      for: props.htmlFor,
      part: props.part,
      tabindex: props.tabIndex,
      readonly: props.readonly ? '' : undefined,
      style: { ...props.style }
    },
    props.children
  );
});

;// CONCATENATED MODULE: ../node_modules/@jupyter/react-components/lib/Select.js



(0,index_js_.provideJupyterDesignSystem)().register((0,index_js_.jpSelect)());

const Select = (0,react_index_js_.forwardRef)((props, forwardedRef) => {
  const ref = (0,react_index_js_.useRef)(null);
  const {
    className,
    autowidth,
    minimal,
    open,
    position,
    autoWidth,
    value,
    displayValue,
    multiple,
    size,
    length,
    options,
    disabled,
    selectedIndex,
    selectedOptions,
    required,
    ...filteredProps
  } = props;

  /** Event listeners - run once */
  useEventListener(ref, 'input', props.onInput);
  useEventListener(ref, 'change', props.onChange);

  /** Properties - run whenever a property has changed */
  useProperties(ref, 'autoWidth', props.autoWidth);
  useProperties(ref, 'value', props.value);
  useProperties(ref, 'displayValue', props.displayValue);
  useProperties(ref, 'multiple', props.multiple);
  useProperties(ref, 'size', props.size);
  useProperties(ref, 'length', props.length);
  useProperties(ref, 'options', props.options);
  useProperties(ref, 'disabled', props.disabled);
  useProperties(ref, 'selectedIndex', props.selectedIndex);
  useProperties(ref, 'selectedOptions', props.selectedOptions);
  useProperties(ref, 'required', props.required);

  /** Methods - uses `useImperativeHandle` hook to pass ref to component */
  (0,react_index_js_.useImperativeHandle)(forwardedRef, () => ref.current, [ref.current]);

  return react_index_js_default().createElement(
    'jp-select',
    {
      ref,
      ...filteredProps,
      position: props.position,
      class: props.className,
      exportparts: props.exportparts,
      for: props.htmlFor,
      part: props.part,
      tabindex: props.tabIndex,
      autowidth: props.autowidth ? '' : undefined,
      minimal: props.minimal ? '' : undefined,
      open: props.open ? '' : undefined,
      style: { ...props.style }
    },
    props.children
  );
});

;// CONCATENATED MODULE: ../node_modules/@jupyter/react-components/lib/Skeleton.js



(0,index_js_.provideJupyterDesignSystem)().register((0,index_js_.jpSkeleton)());

const Skeleton = (0,react_index_js_.forwardRef)((props, forwardedRef) => {
  const ref = (0,react_index_js_.useRef)(null);
  const { className, fill, shape, pattern, shimmer, ...filteredProps } = props;

  /** Properties - run whenever a property has changed */
  useProperties(ref, 'shimmer', props.shimmer);

  /** Methods - uses `useImperativeHandle` hook to pass ref to component */
  (0,react_index_js_.useImperativeHandle)(forwardedRef, () => ref.current, [ref.current]);

  return react_index_js_default().createElement(
    'jp-skeleton',
    {
      ref,
      ...filteredProps,
      fill: props.fill,
      shape: props.shape,
      pattern: props.pattern,
      class: props.className,
      exportparts: props.exportparts,
      for: props.htmlFor,
      part: props.part,
      tabindex: props.tabIndex,
      style: { ...props.style }
    },
    props.children
  );
});

;// CONCATENATED MODULE: ../node_modules/@jupyter/react-components/lib/Slider.js



(0,index_js_.provideJupyterDesignSystem)().register((0,index_js_.jpSlider)());

const Slider = (0,react_index_js_.forwardRef)((props, forwardedRef) => {
  const ref = (0,react_index_js_.useRef)(null);
  const {
    className,
    readonly,
    orientation,
    mode,
    readOnly,
    valueAsNumber,
    valueTextFormatter,
    min,
    max,
    step,
    disabled,
    required,
    ...filteredProps
  } = props;

  /** Event listeners - run once */
  useEventListener(ref, 'change', props.onChange);

  /** Properties - run whenever a property has changed */
  useProperties(ref, 'readOnly', props.readOnly);
  useProperties(ref, 'valueAsNumber', props.valueAsNumber);
  useProperties(ref, 'valueTextFormatter', props.valueTextFormatter);
  useProperties(ref, 'min', props.min);
  useProperties(ref, 'max', props.max);
  useProperties(ref, 'step', props.step);
  useProperties(ref, 'disabled', props.disabled);
  useProperties(ref, 'required', props.required);

  /** Methods - uses `useImperativeHandle` hook to pass ref to component */
  (0,react_index_js_.useImperativeHandle)(forwardedRef, () => ref.current, [ref.current]);

  return react_index_js_default().createElement(
    'jp-slider',
    {
      ref,
      ...filteredProps,
      orientation: props.orientation,
      mode: props.mode,
      class: props.className,
      exportparts: props.exportparts,
      for: props.htmlFor,
      part: props.part,
      tabindex: props.tabIndex,
      readonly: props.readonly ? '' : undefined,
      style: { ...props.style }
    },
    props.children
  );
});

;// CONCATENATED MODULE: ../node_modules/@jupyter/react-components/lib/SliderLabel.js



(0,index_js_.provideJupyterDesignSystem)().register((0,index_js_.jpSliderLabel)());

const SliderLabel = (0,react_index_js_.forwardRef)((props, forwardedRef) => {
  const ref = (0,react_index_js_.useRef)(null);
  const { className, hideMark, disabled, position, ...filteredProps } = props;

  /** Properties - run whenever a property has changed */

  /** Methods - uses `useImperativeHandle` hook to pass ref to component */
  (0,react_index_js_.useImperativeHandle)(forwardedRef, () => ref.current, [ref.current]);

  // Add web component internal classes on top of `className`
  let allClasses = className ?? '';
  if (ref.current?.disabled) {
    allClasses += ' disabled';
  }

  return react_index_js_default().createElement(
    'jp-slider-label',
    {
      ref,
      ...filteredProps,
      position: props.position,
      class: allClasses.trim(),
      exportparts: props.exportparts,
      for: props.htmlFor,
      part: props.part,
      tabindex: props.tabIndex,
      'hide-mark': props.hideMark ? '' : undefined,
      disabled: props.disabled ? '' : undefined,
      style: { ...props.style }
    },
    props.children
  );
});

;// CONCATENATED MODULE: ../node_modules/@jupyter/react-components/lib/Switch.js



(0,index_js_.provideJupyterDesignSystem)().register((0,index_js_.jpSwitch)());

const Switch = (0,react_index_js_.forwardRef)((props, forwardedRef) => {
  const ref = (0,react_index_js_.useRef)(null);
  const {
    className,
    readonly,
    readOnly,
    checked,
    disabled,
    required,
    ...filteredProps
  } = props;

  /** Event listeners - run once */
  useEventListener(ref, 'change', props.onChange);

  /** Properties - run whenever a property has changed */
  useProperties(ref, 'readOnly', props.readOnly);
  useProperties(ref, 'checked', props.checked);
  useProperties(ref, 'disabled', props.disabled);
  useProperties(ref, 'required', props.required);

  /** Methods - uses `useImperativeHandle` hook to pass ref to component */
  (0,react_index_js_.useImperativeHandle)(forwardedRef, () => ref.current, [ref.current]);

  return react_index_js_default().createElement(
    'jp-switch',
    {
      ref,
      ...filteredProps,
      class: props.className,
      exportparts: props.exportparts,
      for: props.htmlFor,
      part: props.part,
      tabindex: props.tabIndex,
      readonly: props.readonly ? '' : undefined,
      style: { ...props.style }
    },
    props.children
  );
});

;// CONCATENATED MODULE: ../node_modules/@jupyter/react-components/lib/Tab.js



(0,index_js_.provideJupyterDesignSystem)().register((0,index_js_.jpTab)());

const Tab = (0,react_index_js_.forwardRef)((props, forwardedRef) => {
  const ref = (0,react_index_js_.useRef)(null);
  const { className, disabled, ...filteredProps } = props;

  /** Properties - run whenever a property has changed */
  useProperties(ref, 'disabled', props.disabled);

  /** Methods - uses `useImperativeHandle` hook to pass ref to component */
  (0,react_index_js_.useImperativeHandle)(forwardedRef, () => ref.current, [ref.current]);

  // Add web component internal classes on top of `className`
  let allClasses = className ?? '';
  if (ref.current?.classList.contains('vertical')) {
    allClasses += ' vertical';
  }

  return react_index_js_default().createElement(
    'jp-tab',
    {
      ref,
      ...filteredProps,
      class: allClasses.trim(),
      exportparts: props.exportparts,
      for: props.htmlFor,
      part: props.part,
      tabindex: props.tabIndex,
      style: { ...props.style }
    },
    props.children
  );
});

;// CONCATENATED MODULE: ../node_modules/@jupyter/react-components/lib/TabPanel.js



(0,index_js_.provideJupyterDesignSystem)().register((0,index_js_.jpTabPanel)());

const TabPanel = (0,react_index_js_.forwardRef)((props, forwardedRef) => {
  const ref = (0,react_index_js_.useRef)(null);
  const { className, ...filteredProps } = props;

  /** Properties - run whenever a property has changed */

  /** Methods - uses `useImperativeHandle` hook to pass ref to component */
  (0,react_index_js_.useImperativeHandle)(forwardedRef, () => ref.current, [ref.current]);

  return react_index_js_default().createElement(
    'jp-tab-panel',
    {
      ref,
      ...filteredProps,
      class: props.className,
      exportparts: props.exportparts,
      for: props.htmlFor,
      part: props.part,
      tabindex: props.tabIndex,
      style: { ...props.style }
    },
    props.children
  );
});

;// CONCATENATED MODULE: ../node_modules/@jupyter/react-components/lib/Tabs.js



(0,index_js_.provideJupyterDesignSystem)().register((0,index_js_.jpTabs)());

const Tabs = (0,react_index_js_.forwardRef)((props, forwardedRef) => {
  const ref = (0,react_index_js_.useRef)(null);
  const {
    className,
    orientation,
    activeid,
    activeindicator,
    activetab,
    ...filteredProps
  } = props;

  /** Event listeners - run once */
  useEventListener(ref, 'change', props.onChange);

  /** Properties - run whenever a property has changed */
  useProperties(ref, 'activeindicator', props.activeindicator);
  useProperties(ref, 'activetab', props.activetab);

  /** Methods - uses `useImperativeHandle` hook to pass ref to component */
  (0,react_index_js_.useImperativeHandle)(forwardedRef, () => ref.current, [ref.current]);

  return react_index_js_default().createElement(
    'jp-tabs',
    {
      ref,
      ...filteredProps,
      orientation: props.orientation,
      activeid: props.activeid,
      class: props.className,
      exportparts: props.exportparts,
      for: props.htmlFor,
      part: props.part,
      tabindex: props.tabIndex,
      style: { ...props.style }
    },
    props.children
  );
});

;// CONCATENATED MODULE: ../node_modules/@jupyter/react-components/lib/TextArea.js



(0,index_js_.provideJupyterDesignSystem)().register((0,index_js_.jpTextArea)());

const TextArea = (0,react_index_js_.forwardRef)((props, forwardedRef) => {
  const ref = (0,react_index_js_.useRef)(null);
  const {
    className,
    appearance,
    resize,
    form,
    list,
    name,
    placeholder,
    readOnly,
    autofocus,
    maxlength,
    minlength,
    cols,
    rows,
    spellcheck,
    disabled,
    required,
    ...filteredProps
  } = props;

  /** Event listeners - run once */
  useEventListener(ref, 'select', props.onSelect);
  useEventListener(ref, 'change', props.onChange);

  /** Properties - run whenever a property has changed */
  useProperties(ref, 'readOnly', props.readOnly);
  useProperties(ref, 'autofocus', props.autofocus);
  useProperties(ref, 'maxlength', props.maxlength);
  useProperties(ref, 'minlength', props.minlength);
  useProperties(ref, 'cols', props.cols);
  useProperties(ref, 'rows', props.rows);
  useProperties(ref, 'spellcheck', props.spellcheck);
  useProperties(ref, 'disabled', props.disabled);
  useProperties(ref, 'required', props.required);

  /** Methods - uses `useImperativeHandle` hook to pass ref to component */
  (0,react_index_js_.useImperativeHandle)(forwardedRef, () => ref.current, [ref.current]);

  return react_index_js_default().createElement(
    'jp-text-area',
    {
      ref,
      ...filteredProps,
      appearance: props.appearance,
      resize: props.resize,
      form: props.form,
      list: props.list,
      name: props.name,
      placeholder: props.placeholder,
      class: props.className,
      exportparts: props.exportparts,
      for: props.htmlFor,
      part: props.part,
      tabindex: props.tabIndex,
      style: { ...props.style }
    },
    props.children
  );
});

;// CONCATENATED MODULE: ../node_modules/@jupyter/react-components/lib/TextField.js



(0,index_js_.provideJupyterDesignSystem)().register((0,index_js_.jpTextField)());

const TextField = (0,react_index_js_.forwardRef)((props, forwardedRef) => {
  const ref = (0,react_index_js_.useRef)(null);
  const {
    className,
    readonly,
    appearance,
    placeholder,
    type,
    list,
    pattern,
    readOnly,
    autofocus,
    maxlength,
    minlength,
    size,
    spellcheck,
    disabled,
    required,
    ...filteredProps
  } = props;

  /** Event listeners - run once */
  useEventListener(ref, 'change', props.onChange);
  useEventListener(ref, 'input', props.onInput);

  /** Properties - run whenever a property has changed */
  useProperties(ref, 'readOnly', props.readOnly);
  useProperties(ref, 'autofocus', props.autofocus);
  useProperties(ref, 'maxlength', props.maxlength);
  useProperties(ref, 'minlength', props.minlength);
  useProperties(ref, 'size', props.size);
  useProperties(ref, 'spellcheck', props.spellcheck);
  useProperties(ref, 'disabled', props.disabled);
  useProperties(ref, 'required', props.required);

  /** Methods - uses `useImperativeHandle` hook to pass ref to component */
  (0,react_index_js_.useImperativeHandle)(forwardedRef, () => ref.current, [ref.current]);

  return react_index_js_default().createElement(
    'jp-text-field',
    {
      ref,
      ...filteredProps,
      appearance: props.appearance,
      placeholder: props.placeholder,
      type: props.type,
      list: props.list,
      pattern: props.pattern,
      class: props.className,
      exportparts: props.exportparts,
      for: props.htmlFor,
      part: props.part,
      tabindex: props.tabIndex,
      readonly: props.readonly ? '' : undefined,
      style: { ...props.style }
    },
    props.children
  );
});

;// CONCATENATED MODULE: ../node_modules/@jupyter/react-components/lib/Toolbar.js



(0,index_js_.provideJupyterDesignSystem)().register((0,index_js_.jpToolbar)());

const Toolbar = (0,react_index_js_.forwardRef)((props, forwardedRef) => {
  const ref = (0,react_index_js_.useRef)(null);
  const { className, ...filteredProps } = props;

  /** Properties - run whenever a property has changed */

  /** Methods - uses `useImperativeHandle` hook to pass ref to component */
  (0,react_index_js_.useImperativeHandle)(forwardedRef, () => ref.current, [ref.current]);

  return react_index_js_default().createElement(
    'jp-toolbar',
    {
      ref,
      ...filteredProps,
      class: props.className,
      exportparts: props.exportparts,
      for: props.htmlFor,
      part: props.part,
      tabindex: props.tabIndex,
      style: { ...props.style }
    },
    props.children
  );
});

;// CONCATENATED MODULE: ../node_modules/@jupyter/react-components/lib/Tooltip.js



(0,index_js_.provideJupyterDesignSystem)().register((0,index_js_.jpTooltip)());

const Tooltip = (0,react_index_js_.forwardRef)((props, forwardedRef) => {
  const ref = (0,react_index_js_.useRef)(null);
  const {
    className,
    horizontalViewportLock,
    verticalViewportLock,
    anchor,
    delay,
    position,
    autoUpdateMode,
    visible,
    anchorElement,
    ...filteredProps
  } = props;

  /** Event listeners - run once */
  useEventListener(ref, 'dismiss', props.onDismiss);

  /** Properties - run whenever a property has changed */
  useProperties(ref, 'visible', props.visible);
  useProperties(ref, 'anchorElement', props.anchorElement);

  /** Methods - uses `useImperativeHandle` hook to pass ref to component */
  (0,react_index_js_.useImperativeHandle)(forwardedRef, () => ref.current, [ref.current]);

  return react_index_js_default().createElement(
    'jp-tooltip',
    {
      ref,
      ...filteredProps,
      anchor: props.anchor,
      delay: props.delay,
      position: props.position,
      'auto-update-mode': props.autoUpdateMode || props['auto-update-mode'],
      class: props.className,
      exportparts: props.exportparts,
      for: props.htmlFor,
      part: props.part,
      tabindex: props.tabIndex,
      'horizontal-viewport-lock': props.horizontalViewportLock ? '' : undefined,
      'vertical-viewport-lock': props.verticalViewportLock ? '' : undefined,
      style: { ...props.style }
    },
    props.children
  );
});

;// CONCATENATED MODULE: ../node_modules/@jupyter/react-components/lib/TreeItem.js



(0,index_js_.provideJupyterDesignSystem)().register((0,index_js_.jpTreeItem)());

const TreeItem = (0,react_index_js_.forwardRef)((props, forwardedRef) => {
  const ref = (0,react_index_js_.useRef)(null);
  const { className, expanded, selected, disabled, ...filteredProps } = props;

  /** Event listeners - run once */
  useEventListener(ref, 'expanded-change', props.onExpand);
  useEventListener(ref, 'selected-change', props.onSelect);

  /** Properties - run whenever a property has changed */
  useProperties(ref, 'expanded', props.expanded);
  useProperties(ref, 'selected', props.selected);
  useProperties(ref, 'disabled', props.disabled);

  /** Methods - uses `useImperativeHandle` hook to pass ref to component */
  (0,react_index_js_.useImperativeHandle)(forwardedRef, () => ref.current, [ref.current]);

  // Add web component internal classes on top of `className`
  let allClasses = className ?? '';
  if (ref.current?.nested) {
    allClasses += ' nested';
  }

  return react_index_js_default().createElement(
    'jp-tree-item',
    {
      ref,
      ...filteredProps,
      class: allClasses.trim(),
      exportparts: props.exportparts,
      for: props.htmlFor,
      part: props.part,
      tabindex: props.tabIndex,
      style: { ...props.style }
    },
    props.children
  );
});

;// CONCATENATED MODULE: ../node_modules/@jupyter/react-components/lib/TreeView.js




(0,index_js_.provideJupyterDesignSystem)().register((0,index_js_.jpTreeView)());

const TreeView = (0,react_index_js_.forwardRef)((props, forwardedRef) => {
  const ref = (0,react_index_js_.useRef)(null);
  const { className, renderCollapsedNodes, currentSelected, ...filteredProps } =
    props;

  (0,react_index_js_.useLayoutEffect)(() => {
    // Fix using private API to force refresh of nested flag on
    // first level of tree items.
    ref.current?.setItems();
  }, [ref.current]);

  /** Properties - run whenever a property has changed */
  useProperties(ref, 'currentSelected', props.currentSelected);

  /** Methods - uses `useImperativeHandle` hook to pass ref to component */
  (0,react_index_js_.useImperativeHandle)(forwardedRef, () => ref.current, [ref.current]);

  return react_index_js_default().createElement(
    'jp-tree-view',
    {
      ref,
      ...filteredProps,
      class: props.className,
      exportparts: props.exportparts,
      for: props.htmlFor,
      part: props.part,
      tabindex: props.tabIndex,
      'render-collapsed-nodes': props.renderCollapsedNodes ? '' : undefined,
      style: { ...props.style }
    },
    props.children
  );
});

;// CONCATENATED MODULE: ../node_modules/@jupyter/react-components/lib/Picker.js



(0,index_js_.provideJupyterDesignSystem)().register((0,index_js_.jpPicker)());

const Picker = (0,react_index_js_.forwardRef)((props, forwardedRef) => {
  const ref = (0,react_index_js_.useRef)(null);
  const {
    className,
    filterSelected,
    filterQuery,
    selection,
    options,
    maxSelected,
    noSuggestionsText,
    suggestionsAvailableText,
    loadingText,
    label,
    labelledby,
    placeholder,
    menuPlacement,
    showLoading,
    listItemTemplate,
    defaultListItemTemplate,
    menuOptionTemplate,
    defaultMenuOptionTemplate,
    listItemContentsTemplate,
    menuOptionContentsTemplate,
    optionsList,
    query,
    itemsPlaceholderElement,
    ...filteredProps
  } = props;

  /** Properties - run whenever a property has changed */
  useProperties(ref, 'showLoading', props.showLoading);
  useProperties(ref, 'listItemTemplate', props.listItemTemplate);
  useProperties(ref, 'defaultListItemTemplate', props.defaultListItemTemplate);
  useProperties(ref, 'menuOptionTemplate', props.menuOptionTemplate);
  useProperties(
    ref,
    'defaultMenuOptionTemplate',
    props.defaultMenuOptionTemplate
  );
  useProperties(
    ref,
    'listItemContentsTemplate',
    props.listItemContentsTemplate
  );
  useProperties(
    ref,
    'menuOptionContentsTemplate',
    props.menuOptionContentsTemplate
  );
  useProperties(ref, 'optionsList', props.optionsList);
  useProperties(ref, 'query', props.query);
  useProperties(ref, 'itemsPlaceholderElement', props.itemsPlaceholderElement);

  /** Methods - uses `useImperativeHandle` hook to pass ref to component */
  (0,react_index_js_.useImperativeHandle)(forwardedRef, () => ref.current, [ref.current]);

  return react_index_js_default().createElement(
    'jp-draft-picker',
    {
      ref,
      ...filteredProps,
      selection: props.selection,
      options: props.options,
      'max-selected': props.maxSelected || props['max-selected'],
      'no-suggestions-text':
        props.noSuggestionsText || props['no-suggestions-text'],
      'suggestions-available-text':
        props.suggestionsAvailableText || props['suggestions-available-text'],
      'loading-text': props.loadingText || props['loading-text'],
      label: props.label,
      labelledby: props.labelledby,
      placeholder: props.placeholder,
      'menu-placement': props.menuPlacement || props['menu-placement'],
      class: props.className,
      exportparts: props.exportparts,
      for: props.htmlFor,
      part: props.part,
      tabindex: props.tabIndex,
      'filter-selected': props.filterSelected ? '' : undefined,
      'filter-query': props.filterQuery ? '' : undefined,
      style: { ...props.style }
    },
    props.children
  );
});

;// CONCATENATED MODULE: ../node_modules/@jupyter/react-components/lib/PickerMenu.js



(0,index_js_.provideJupyterDesignSystem)().register((0,index_js_.jpPickerMenu)());

const PickerMenu = (0,react_index_js_.forwardRef)((props, forwardedRef) => {
  const ref = (0,react_index_js_.useRef)(null);
  const { className, suggestionsAvailableText, ...filteredProps } = props;

  /** Properties - run whenever a property has changed */
  useProperties(
    ref,
    'suggestionsAvailableText',
    props.suggestionsAvailableText
  );

  /** Methods - uses `useImperativeHandle` hook to pass ref to component */
  (0,react_index_js_.useImperativeHandle)(forwardedRef, () => ref.current, [ref.current]);

  return react_index_js_default().createElement(
    'jp-draft-picker-menu',
    {
      ref,
      ...filteredProps,
      class: props.className,
      exportparts: props.exportparts,
      for: props.htmlFor,
      part: props.part,
      tabindex: props.tabIndex,
      style: { ...props.style }
    },
    props.children
  );
});

;// CONCATENATED MODULE: ../node_modules/@jupyter/react-components/lib/PickerMenuOption.js



(0,index_js_.provideJupyterDesignSystem)().register((0,index_js_.jpPickerMenuOption)());

const PickerMenuOption = (0,react_index_js_.forwardRef)((props, forwardedRef) => {
  const ref = (0,react_index_js_.useRef)(null);
  const { className, value, contentsTemplate, ...filteredProps } = props;

  /** Properties - run whenever a property has changed */
  useProperties(ref, 'contentsTemplate', props.contentsTemplate);

  /** Methods - uses `useImperativeHandle` hook to pass ref to component */
  (0,react_index_js_.useImperativeHandle)(forwardedRef, () => ref.current, [ref.current]);

  return react_index_js_default().createElement(
    'jp-draft-picker-menu-option',
    {
      ref,
      ...filteredProps,
      value: props.value,
      class: props.className,
      exportparts: props.exportparts,
      for: props.htmlFor,
      part: props.part,
      tabindex: props.tabIndex,
      style: { ...props.style }
    },
    props.children
  );
});

;// CONCATENATED MODULE: ../node_modules/@jupyter/react-components/lib/PickerList.js



(0,index_js_.provideJupyterDesignSystem)().register((0,index_js_.jpPickerList)());

const PickerList = (0,react_index_js_.forwardRef)((props, forwardedRef) => {
  const ref = (0,react_index_js_.useRef)(null);
  const { className, ...filteredProps } = props;

  /** Properties - run whenever a property has changed */

  /** Methods - uses `useImperativeHandle` hook to pass ref to component */
  (0,react_index_js_.useImperativeHandle)(forwardedRef, () => ref.current, [ref.current]);

  return react_index_js_default().createElement(
    'jp-draft-picker-list',
    {
      ref,
      ...filteredProps,
      class: props.className,
      exportparts: props.exportparts,
      for: props.htmlFor,
      part: props.part,
      tabindex: props.tabIndex,
      style: { ...props.style }
    },
    props.children
  );
});

;// CONCATENATED MODULE: ../node_modules/@jupyter/react-components/lib/PickerListItem.js



(0,index_js_.provideJupyterDesignSystem)().register((0,index_js_.jpPickerListItem)());

const PickerListItem = (0,react_index_js_.forwardRef)((props, forwardedRef) => {
  const ref = (0,react_index_js_.useRef)(null);
  const { className, value, contentsTemplate, ...filteredProps } = props;

  /** Properties - run whenever a property has changed */
  useProperties(ref, 'contentsTemplate', props.contentsTemplate);

  /** Methods - uses `useImperativeHandle` hook to pass ref to component */
  (0,react_index_js_.useImperativeHandle)(forwardedRef, () => ref.current, [ref.current]);

  return react_index_js_default().createElement(
    'jp-draft-picker-list-item',
    {
      ref,
      ...filteredProps,
      value: props.value,
      class: props.className,
      exportparts: props.exportparts,
      for: props.htmlFor,
      part: props.part,
      tabindex: props.tabIndex,
      style: { ...props.style }
    },
    props.children
  );
});

;// CONCATENATED MODULE: ../node_modules/@jupyter/react-components/lib/index.js


















































/***/ })

}]);
//# sourceMappingURL=2816.03541f3103bf4c09e591.js.map?v=03541f3103bf4c09e591